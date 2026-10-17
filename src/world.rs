//! The world: grid, molecules, binding, the timestep (spec: The world, Binding, One timestep).

use crate::align::{bind_prob, Aligner};
use crate::config::Config;
use crate::reaction::{exec_step, ChemParams, Ctx, Outcome, Pointers, Strand};
use crate::rng::{Rng, Xoshiro};
use crate::symbols::{complement, Table};
use std::collections::HashMap;

/// Event counters, reset by whoever reads them.
#[derive(Clone, Debug, Default)]
pub struct Stats {
    pub encounters: u64,
    pub binds: u64,
    pub kin_binds: u64,
    pub ended: u64,
    pub products: u64,
    pub pair_decays: u64,
    pub displacements: u64,
    pub instr: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Bond {
    Free,
    /// The executing molecule: its partner, pointers, and both species at binding time.
    /// p0: binding probability when the bond formed; copied: `=` instructions run since.
    Active { partner: usize, ptr: Pointers, sp_act: u32, sp_pas: u32, p0: f32, copied: u32 },
    Passive { partner: usize },
}

#[derive(Clone, Debug)]
pub struct Mol {
    pub strand: Strand,
    pub x: usize,
    pub y: usize,
    pub species: u32,
    pub bond: Bond,
    pub alive: bool,
}

/// What happened during a step that the output layer needs to know about.
pub trait Observer {
    fn new_species(&mut self, id: u32, born: u64, pa: i64, pp: i64, seq: &[u8]);
    fn product(&mut self, sp_act: u32, sp_pas: u32, sp_prod: u32);
}

pub struct World {
    pub gx: usize,
    pub gy: usize,
    /// cell -> molecule index + 1, 0 for empty
    pub grid: Vec<u32>,
    pub mols: Vec<Mol>,
    pub free_slots: Vec<usize>,
    pub species: HashMap<Vec<u8>, u32>,
    pub next_species: u32,
    pub step: u64,
    pub rng: Xoshiro,
    pub decay: f64,
    pub bindlen: usize,
    pub compat_partner: bool,
    pub decay_independent: bool,
    pub encounter_turn: bool,
    pub displace: bool,
    pub displace_k: f64,
    pub displace_rate: f64,
    pub chem: ChemParams,
    pub stats: Stats,
    tab: Table,
    al: Aligner,
    acted: Vec<u64>,
    order: Vec<usize>,
    comp: Vec<u8>,
    products: Vec<Vec<u8>>,
    label: Vec<u8>,
    other: Vec<u8>,
}

impl World {
    pub fn new(cfg: &Config, obs: &mut dyn Observer) -> Result<World, String> {
        let mut w = World {
            gx: cfg.gridx,
            gy: cfg.gridy,
            grid: vec![0; cfg.gridx * cfg.gridy],
            mols: Vec::new(),
            free_slots: Vec::new(),
            species: HashMap::new(),
            next_species: 1,
            step: 0,
            rng: Xoshiro::new(cfg.seed),
            decay: cfg.decay,
            bindlen: cfg.bindlen,
            compat_partner: cfg.compat_partner,
            decay_independent: cfg.decay_independent,
            encounter_turn: cfg.encounter_turn,
            displace: cfg.displace,
            displace_k: cfg.displace_k,
            displace_rate: cfg.displace_rate,
            chem: ChemParams { maxl: cfg.maxl, mutate: cfg.mutate, compat_deletion_skip: false },
            stats: Stats::default(),
            tab: Table::new(),
            al: Aligner::new(),
            acted: Vec::new(),
            order: Vec::new(),
            comp: Vec::new(),
            products: Vec::new(),
            label: Vec::new(),
            other: Vec::new(),
        };
        for (seq, x, y) in &cfg.agents {
            if w.grid[w.cell(*x, *y)] != 0 {
                return Err(format!("two molecules at {x} {y}"));
            }
            let sp = w.species_of(seq, -1, -1, obs);
            w.add(Strand::new(seq, cfg.maxl), *x, *y, sp);
        }
        Ok(w)
    }

    #[inline]
    pub fn cell(&self, x: usize, y: usize) -> usize {
        y * self.gx + x
    }

    pub fn count(&self) -> usize {
        self.mols.len() - self.free_slots.len()
    }

    fn add(&mut self, strand: Strand, x: usize, y: usize, species: u32) -> usize {
        let m = Mol { strand, x, y, species, bond: Bond::Free, alive: true };
        let idx = if let Some(i) = self.free_slots.pop() {
            self.mols[i] = m;
            i
        } else {
            self.mols.push(m);
            self.acted.push(0);
            self.mols.len() - 1
        };
        let c = self.cell(x, y);
        self.grid[c] = idx as u32 + 1;
        self.acted[idx] = self.step + 1; // products don't act in the step they're made
        idx
    }

    fn remove(&mut self, idx: usize) {
        let c = self.cell(self.mols[idx].x, self.mols[idx].y);
        self.grid[c] = 0;
        self.mols[idx].alive = false;
        self.mols[idx].bond = Bond::Free;
        self.free_slots.push(idx);
    }

    /// Species id for a sequence, creating it (and reporting it) if new.
    fn species_of(&mut self, seq: &[u8], pa: i64, pp: i64, obs: &mut dyn Observer) -> u32 {
        if let Some(&id) = self.species.get(seq) {
            return id;
        }
        let id = self.next_species;
        self.next_species += 1;
        self.species.insert(seq.to_vec(), id);
        obs.new_species(id, self.step + 1, pa, pp, seq);
        id
    }

    /// Restore an entry from species.tsv when resuming.
    pub fn known_species(&mut self, id: u32, seq: &[u8]) {
        self.species.insert(seq.to_vec(), id);
        if id >= self.next_species {
            self.next_species = id + 1;
        }
    }

    pub fn neighbours(&self, x: usize, y: usize) -> [(usize, usize); 8] {
        let (gx, gy) = (self.gx, self.gy);
        let mut out = [(0, 0); 8];
        let mut k = 0;
        for dy in [gy - 1, 0, 1] {
            for dx in [gx - 1, 0, 1] {
                if dx == 0 && dy == 0 {
                    continue;
                }
                out[k] = ((x + dx) % gx, (y + dy) % gy);
                k += 1;
            }
        }
        out
    }

    /// The old code's neighbour scan order: x offset outer, y offset inner, -1..=1.
    fn scan_order(&self, x: usize, y: usize) -> [(usize, usize); 8] {
        let (gx, gy) = (self.gx, self.gy);
        let mut out = [(0, 0); 8];
        let mut k = 0;
        for dx in [gx - 1, 0, 1] {
            for dy in [gy - 1, 0, 1] {
                if dx == 0 && dy == 0 {
                    continue;
                }
                out[k] = ((x + dx) % gx, (y + dy) % gy);
                k += 1;
            }
        }
        out
    }

    /// End a molecule's reaction from its side: it becomes free, and its species is
    /// re-read from its (possibly changed) sequence.
    fn release(&mut self, idx: usize, sp_act: u32, sp_pas: u32, obs: &mut dyn Observer) {
        self.mols[idx].bond = Bond::Free;
        let seq = self.mols[idx].strand.seq().to_vec();
        let sp = self.species_of(&seq, sp_act as i64, sp_pas as i64, obs);
        self.mols[idx].species = sp;
    }

    fn place_product(&mut self, seq: Vec<u8>, near: usize, sp_act: u32, sp_pas: u32, obs: &mut dyn Observer) {
        let (x, y) = (self.mols[near].x, self.mols[near].y);
        let mut empty = [(0usize, 0usize); 8];
        let mut n = 0;
        for (nx, ny) in self.neighbours(x, y) {
            if self.grid[self.cell(nx, ny)] == 0 {
                empty[n] = (nx, ny);
                n += 1;
            }
        }
        let sp = self.species_of(&seq, sp_act as i64, sp_pas as i64, obs);
        obs.product(sp_act, sp_pas, sp);
        if n == 0 {
            return; // no room: the product is lost
        }
        let (px, py) = empty[self.rng.below(n)];
        let maxl = self.chem.maxl;
        self.add(Strand::new(&seq, maxl), px, py, sp);
    }

    fn pair_of(&self, idx: usize) -> Option<(usize, usize)> {
        match self.mols[idx].bond {
            Bond::Free => None,
            Bond::Active { partner, .. } => Some((idx, partner)),
            Bond::Passive { partner } => Some((partner, idx)),
        }
    }

    pub fn step(&mut self, obs: &mut dyn Observer) {
        let now = self.step + 1;
        self.order.clear();
        for (i, m) in self.mols.iter().enumerate() {
            if m.alive {
                self.order.push(i);
            }
        }
        // Fisher-Yates
        for i in (1..self.order.len()).rev() {
            let j = self.rng.below(i + 1);
            self.order.swap(i, j);
        }
        let order = std::mem::take(&mut self.order);
        for &idx in &order {
            if !self.mols[idx].alive || self.acted[idx] >= now {
                continue;
            }
            self.acted[idx] = now;
            match self.pair_of(idx) {
                Some((a, p)) => self.act_pair(idx, a, p, now, obs),
                None => self.act_free(idx, now, obs),
            }
        }
        self.order = order;
        self.step = now;
    }

    fn roll_decay(&mut self) -> bool {
        self.rng.next_f64() < self.decay
    }

    fn act_free(&mut self, idx: usize, now: u64, obs: &mut dyn Observer) {
        if self.roll_decay() {
            self.remove(idx);
            return;
        }
        let (x, y) = (self.mols[idx].x, self.mols[idx].y);
        let mut cand = [0usize; 8];
        let mut n = 0;
        let cells = if self.compat_partner { self.scan_order(x, y) } else { self.neighbours(x, y) };
        for (nx, ny) in cells {
            let g = self.grid[self.cell(nx, ny)];
            if g != 0 {
                let j = (g - 1) as usize;
                let free = self.mols[j].bond == Bond::Free;
                if (free && self.acted[j] < now) || (!free && self.displace) {
                    cand[n] = j;
                    n += 1;
                }
            }
        }
        if n == 0 {
            return;
        }
        let b = if self.compat_partner {
            self.rng.next_f64(); // the old code drew a number it then ignored
            cand[0]
        } else {
            cand[self.rng.below(n)]
        };
        let busy = self.mols[b].bond != Bond::Free;
        if self.encounter_turn && !busy {
            self.acted[b] = now; // an encounter uses up the partner's turn, bound or not
        }
        self.stats.encounters += 1;
        self.comp.clear();
        let a_seq = self.mols[idx].strand.seq();
        self.comp.extend(a_seq.iter().map(|&c| complement(c)));
        let sw = self.al.align(&self.comp, self.mols[b].strand.seq(), &self.tab);
        let p = bind_prob(&sw, self.bindlen);
        let rno = self.rng.next_f64() as f32;
        if busy {
            // strand displacement: the newcomer must bind, and the old bond must give way.
            // The old bond starts at its formation probability p0 and hardens as its copy
            // grows; after DISPLACE_K copy instructions it can no longer be broken into.
            let (oa, op) = self.pair_of(b).unwrap();
            let (p0, copied, sa, sp) = match self.mols[oa].bond {
                Bond::Active { p0, copied, sp_act, sp_pas, .. } => (p0, copied, sp_act, sp_pas),
                _ => unreachable!(),
            };
            let soft = (1.0 - copied as f64 / self.displace_k).max(0.0);
            let take = self.displace_rate * p as f64 * (1.0 - p0 as f64) * soft;
            if (rno as f64) >= take {
                return;
            }
            let other = if b == oa { op } else { oa };
            self.release(oa, sa, sp, obs);
            self.release(op, sa, sp, obs);
            self.acted[other] = now;
            self.stats.displacements += 1;
        } else if rno >= p {
            return;
        }
        let (act, pas, act_idx, pas_idx) = if sw.s1 >= sw.s2 { (idx, b, sw.s1, sw.s2) } else { (b, idx, sw.s2, sw.s1) };
        let (sp_act, sp_pas) = (self.mols[act].species, self.mols[pas].species);
        self.stats.binds += 1;
        self.acted[b] = now;
        if sp_act == sp_pas {
            self.stats.kin_binds += 1;
        }
        self.mols[act].bond = Bond::Active { partner: pas, ptr: Pointers::at_bind(act_idx, pas_idx), sp_act, sp_pas, p0: p, copied: 0 };
        self.mols[pas].bond = Bond::Passive { partner: act };
    }

    fn act_pair(&mut self, drawn: usize, act: usize, pas: usize, now: u64, obs: &mut dyn Observer) {
        let other = if drawn == act { pas } else { act };
        self.acted[other] = now;
        let (sp_act, sp_pas) = match self.mols[act].bond {
            Bond::Active { sp_act, sp_pas, .. } => (sp_act, sp_pas),
            _ => unreachable!(),
        };
        if !self.decay_independent {
            // baseline (old) rule: one roll for the pair, both die
            if self.roll_decay() {
                self.stats.pair_decays += 1;
                self.remove(drawn);
                self.remove(other);
                return;
            }
        } else {
        // DECAY_INDEPENDENT: each member rolls; a decaying molecule dies, its partner is freed
        let d_drawn = self.roll_decay();
        let d_other = self.roll_decay();
        if d_drawn || d_other {
            for (m, dies) in [(drawn, d_drawn), (other, d_other)] {
                if dies {
                    self.remove(m);
                } else {
                    self.release(m, sp_act, sp_pas, obs);
                }
            }
            return;
        }
        }
        // one instruction of the active's program
        let mut ptr = match &self.mols[act].bond {
            Bond::Active { ptr, .. } => ptr.clone(),
            _ => unreachable!(),
        };
        if self.displace {
            let at = ptr.pos[crate::reaction::I][ptr.side[crate::reaction::I]];
            let c = if ptr.side[crate::reaction::I] == 1 { self.mols[act].strand.at(at) } else { self.mols[pas].strand.at(at) };
            if c == b'=' {
                if let Bond::Active { copied, .. } = &mut self.mols[act].bond {
                    *copied += 1;
                }
            }
        }
        let (own, partner) = two_mut(&mut self.mols, act, pas);
        let mut ctx = Ctx {
            tab: &self.tab,
            al: &mut self.al,
            params: self.chem,
            label: std::mem::take(&mut self.label),
            other: std::mem::take(&mut self.other),
        };
        self.products.clear();
        let out = exec_step(&mut ptr, &mut own.strand, &mut partner.strand, &mut ctx, &mut self.rng, &mut self.products);
        self.label = ctx.label;
        self.other = ctx.other;
        if let Bond::Active { ptr: p, .. } = &mut self.mols[act].bond {
            *p = ptr;
        }
        self.stats.instr += 1;
        self.stats.products += self.products.len() as u64;
        let products = std::mem::take(&mut self.products);
        for seq in products.iter() {
            self.place_product(seq.clone(), act, sp_act, sp_pas, obs);
        }
        self.products = products;
        match out {
            Outcome::Continue => {}
            Outcome::Ended => {
                self.stats.ended += 1;
                self.release(act, sp_act, sp_pas, obs);
                self.release(pas, sp_act, sp_pas, obs);
            }
            Outcome::ActiveDestroyed => {
                self.remove(act);
                self.release(pas, sp_act, sp_pas, obs);
            }
            Outcome::PassiveDestroyed => {
                self.remove(pas);
                self.release(act, sp_act, sp_pas, obs);
            }
        }
    }
}

fn two_mut(v: &mut [Mol], a: usize, b: usize) -> (&mut Mol, &mut Mol) {
    assert_ne!(a, b);
    if a < b {
        let (l, r) = v.split_at_mut(b);
        (&mut l[a], &mut r[0])
    } else {
        let (l, r) = v.split_at_mut(a);
        (&mut r[0], &mut l[b])
    }
}

impl World {
    /// Rebuild a world from a checkpoint written by `output::write_checkpoint`.
    /// `cfg` supplies the parameters (the checkpoint's own CONFIG block, normally).
    pub fn restore(cfg: &Config, text: &str) -> Result<World, String> {
        struct Nop;
        impl Observer for Nop {
            fn new_species(&mut self, _: u32, _: u64, _: i64, _: i64, _: &[u8]) {}
            fn product(&mut self, _: u32, _: u32, _: u32) {}
        }
        let mut empty = cfg.clone();
        empty.agents.clear();
        let mut w = World::new(&empty, &mut Nop)?;
        let mut lines = text.lines();
        let mut next = || lines.next().ok_or("truncated checkpoint".to_string());
        let hdr = next()?;
        if hdr != "SMOLCKPT 1" {
            return Err(format!("not a checkpoint: {hdr}"));
        }
        let field = |l: &str, k: &str| -> Result<Vec<u64>, String> {
            let mut it = l.split_whitespace();
            if it.next() != Some(k) {
                return Err(format!("expected {k}, got {l}"));
            }
            it.map(|v| v.parse::<u64>().map_err(|e| e.to_string())).collect()
        };
        w.step = field(next()?, "STEP")?[0];
        let r = field(next()?, "RNG")?;
        w.rng.s = [r[0], r[1], r[2], r[3]];
        w.next_species = field(next()?, "NEXTSPECIES")?[0] as u32;
        let nconf = field(next()?, "CONFIG")?[0];
        for _ in 0..nconf {
            next()?;
        }
        let slots = field(next()?, "SLOTS")?[0] as usize;
        let free = field(next()?, "FREE")?;
        let nm = field(next()?, "MOLS")?[0] as usize;
        let maxl = cfg.maxl;
        for _ in 0..nm {
            let l = next()?;
            let p: Vec<&str> = l.split_whitespace().collect();
            let idx: usize = p[0].parse().map_err(|_| "bad index")?;
            let (x, y): (usize, usize) = (p[1].parse().unwrap(), p[2].parse().unwrap());
            let species: u32 = p[3].parse().unwrap();
            let seq: Vec<u8> = p[4].bytes().map(|c| if c == b'.' { 0 } else { c }).collect();
            let mut strand = Strand::new(&[], maxl);
            strand.buf[..seq.len()].copy_from_slice(&seq);
            let bond = match p[5] {
                "F" => Bond::Free,
                "P" => Bond::Passive { partner: p[6].parse().unwrap() },
                "A" => {
                    let n: Vec<usize> = p[6..21].iter().map(|v| v.parse().unwrap()).collect();
                    let mut ptr = Pointers::default();
                    for k in 0..4 {
                        ptr.pos[k] = [n[3 + 2 * k], n[4 + 2 * k]];
                        ptr.side[k] = n[11 + k];
                    }
                    let (p0, copied) = if p.len() > 21 { (p[21].parse().unwrap(), p[22].parse().unwrap()) } else { (0.0, 0) };
                    Bond::Active { partner: n[0], ptr, sp_act: n[1] as u32, sp_pas: n[2] as u32, p0, copied }
                }
                b => return Err(format!("bad bond {b}")),
            };
            while w.mols.len() <= idx {
                w.mols.push(Mol { strand: Strand::new(&[], 0), x: 0, y: 0, species: 0, bond: Bond::Free, alive: false });
                w.acted.push(0);
            }
            let c = w.cell(x, y);
            w.grid[c] = idx as u32 + 1;
            w.mols[idx] = Mol { strand, x, y, species, bond, alive: true };
        }
        while w.mols.len() < slots {
            w.mols.push(Mol { strand: Strand::new(&[], 0), x: 0, y: 0, species: 0, bond: Bond::Free, alive: false });
            w.acted.push(0);
        }
        // free list in its saved order, so a resumed run continues exactly
        w.free_slots = free[1..].iter().map(|&v| v as usize).collect();
        Ok(w)
    }
}
