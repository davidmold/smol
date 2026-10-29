//! Config files, in the old Stringmol format: `KEY value` lines, `%` or `#` comments,
//! `AGENT seq [count] ...` followed by `GRIDPOS x y` lines.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Config {
    pub gridx: usize,
    pub gridy: usize,
    pub decay: f64,
    /// DAMAGE: chance per symbol per step that a symbol is replaced by a random one, in any
    /// molecule, free or bound (an environmental hazard, like UV).
    pub damage: f64,
    /// DAMAGE_SHIELD_BOUND: molecules in a reaction are shielded from damage.
    pub damage_shield_bound: bool,
    /// RATE_MIN / RATE_MAX: reaction-rate field across x. A molecule's turn (binding attempt or
    /// one instruction) happens with probability rate(x); decay is unaffected. The field runs
    /// from RATE_MIN at x = 0 to RATE_MAX at the middle column and back (smooth on the torus).
    pub rate_min: f64,
    pub rate_max: f64,
    /// VERTICAL: a side-on world. y = gy-1 is the ground, y = 0 the top; no wrap in y. A product
    /// may only be placed on the ground or on top of an occupied cell (support); when a molecule
    /// is removed, the column above it falls by one. RATE_MIN /
    /// RATE_MAX run from the ground up (light) instead of across x.
    pub vertical: bool,
    /// DEPTH n > 1: n slices of GRIDX x GRIDY side by side (a thin 3D slab; no wrap between
    /// slices). GRIDPOS takes an optional slice number z. Neighbours are the 3x3x3 block.
    pub depth: usize,
    /// SHADE s > 0 (VERTICAL): a molecule's reaction rate is multiplied by s for every molecule
    /// stacked above it in its column (on top of any RATE_MIN/RATE_MAX light-by-height field).
    pub shade: f64,
    /// DECAY_CHEM S > 0: each symbol gets a stability weight S^u (u uniform in -1..1, table from
    /// DECAY_CHEM_SEED), scaled so the alphabet's mean weight is 1; a molecule's decay chance is
    /// DECAY times the mean weight of its symbols.
    pub decay_chem: f64,
    pub decay_chem_seed: u64,
    pub mutate: f32,
    pub maxl: usize,
    pub seed: u64,
    pub nsteps: u64,
    pub popdy_every: u64,
    pub snap_every: u64,
    pub copy_every: u64,
    pub ckpt_every: u64,
    pub keep_ckpt: usize,
    /// LIVEEVERY: write live.grid and live.json (latest state, overwritten) every N steps; 0 = off
    pub live_every: u64,
    pub bindlen: usize,
    /// > 0: absolute-score binding rule with this scale (see align::bind_prob_rule)
    pub bind_t: f64,
    /// Old partner choice: always the first free neighbour in scan order (quirk 1).
    pub compat_partner: bool,
    /// Decay in pairs. Default (as the old code): one roll per pair, both die.
    /// DECAY_INDEPENDENT 1: each member rolls; a decaying molecule dies, its partner is freed.
    pub decay_independent: bool,
    /// A failed encounter also uses up the partner's turn (old rule; default 0 = it does not).
    pub encounter_turn: bool,
    /// Strand displacement (design note): a free molecule may break into a reaction.
    pub displace: bool,
    /// Copy instructions after which a reaction can no longer be broken into.
    pub displace_k: f64,
    /// Overall scale on the displacement probability.
    pub displace_rate: f64,
    /// On displacement, undo the abandoned reaction's writes to both strands.
    pub displace_restore: bool,
    /// Two-sided reactions: the passive also runs its own code, one instruction after the
    /// active's each step.
    pub twosided: bool,
    /// Copies are built on a separate nascent strand (see reaction::copy_nascent).
    pub strandcopy: bool,
    /// OFFSPRINGLOG: copies.tsv gains a column of placed products (the count column is all
    /// attempts, as before), and lives.tsv records, per species and window, the molecules that
    /// died: how long they lived, how many placed products they made, and how many placed exact
    /// copies of them were made (their successful offspring).
    pub offspring_log: bool,
    /// SPECIESFORGET 1: at each checkpoint forget species absent from the world (memory)
    pub species_forget: bool,
    /// DISTURB_*: every DISTURB_EVERY steps (0 = off) up to DISTURB_UNTIL (0 = no end), remove
    /// every molecule in DISTURB_SIZE^2 cells: one random square patch (DISTURB_MODE 1) or that
    /// many random cells across the world (DISTURB_MODE 2). Cells come from DISTURB_SEED and the step.
    pub disturb_every: u64,
    pub disturb_until: u64,
    pub disturb_size: usize,
    pub disturb_mode: u32,
    pub disturb_seed: u64,
    /// GZIPLOGS: popdy001.dat, species.tsv, copies.tsv and lives.tsv are written gzipped (.gz),
    /// one gzip member per checkpoint interval. Doesn't change the run.
    pub gziplogs: bool,
    /// CKPTKEEPEVERY n > 0: checkpoints at multiples of n are kept permanently.
    pub ckpt_keep_every: u64,
    /// TRANSLATE: `^C` switches `=` to translating letter pairs through a codon table.
    pub translate: bool,
    pub translate_seed: u64,
    pub translate_stop: bool,
    /// (sequence, x, y)
    pub agents: Vec<(Vec<u8>, usize, usize)>,
    /// Every key/value as read, for run.json.
    pub raw: Vec<(String, String)>,
}

fn strip(line: &str) -> &str {
    let l = line.trim();
    if l.starts_with('%') || l.starts_with('#') {
        ""
    } else {
        l
    }
}

impl Config {
    pub fn parse(text: &str) -> Result<Config, String> {
        let mut kv: HashMap<String, String> = HashMap::new();
        let mut raw = Vec::new();
        let mut agents = Vec::new();
        let mut agent_z: Vec<usize> = Vec::new();
        let mut cur: Option<Vec<u8>> = None;
        for line in text.lines() {
            let l = strip(line);
            if l.is_empty() {
                continue;
            }
            let mut it = l.split_whitespace();
            let key = it.next().unwrap().to_string();
            let rest: Vec<&str> = it.collect();
            match key.as_str() {
                "AGENT" => {
                    let seq = rest.first().ok_or("AGENT without a sequence")?;
                    cur = Some(seq.as_bytes().to_vec());
                }
                "GRIDPOS" => {
                    let seq = cur.clone().ok_or("GRIDPOS before any AGENT")?;
                    let x: usize = rest.first().and_then(|v| v.parse().ok()).ok_or("bad GRIDPOS")?;
                    let y: usize = rest.get(1).and_then(|v| v.parse().ok()).ok_or("bad GRIDPOS")?;
                    let z: usize = rest.get(2).and_then(|v| v.parse().ok()).unwrap_or(0);
                    agents.push((seq, x, y));
                    agent_z.push(z);
                }
                _ => {
                    if let Some(v) = rest.first() {
                        kv.insert(key.clone(), v.to_string());
                        raw.push((key, v.to_string()));
                    }
                }
            }
        }
        let get = |k: &str| kv.get(k).cloned();
        let num = |k: &str, d: f64| -> Result<f64, String> {
            match get(k) {
                Some(v) => v.parse::<f64>().map_err(|_| format!("bad value for {k}: {v}")),
                None => Ok(d),
            }
        };
        let mutate = match get("MUTATE") {
            Some(v) => v.parse::<f32>().map_err(|_| format!("bad MUTATE {v}"))?,
            None => return Err("MUTATE is required (spec: Copying and mutation)".into()),
        };
        let gridx = num("GRIDX", 125.0)? as usize;
        let gridy = num("GRIDY", 100.0)? as usize;
        let depth = (num("DEPTH", 1.0)? as usize).max(1);
        for (i, (_, x, y)) in agents.iter_mut().enumerate() {
            if *x >= gridx || *y >= gridy || agent_z[i] >= depth {
                return Err(format!("GRIDPOS {x} {y} {} is outside the {gridx}x{gridy}x{depth} grid", agent_z[i]));
            }
            *x += agent_z[i] * gridx; // slices sit side by side in storage
        }
        Ok(Config {
            gridx,
            gridy,
            decay: num("DECAY", 0.0005)?,
            damage: num("DAMAGE", 0.0)?,
            damage_shield_bound: num("DAMAGE_SHIELD_BOUND", 0.0)? != 0.0,
            rate_min: num("RATE_MIN", 1.0)?,
            rate_max: num("RATE_MAX", 1.0)?,
            vertical: num("VERTICAL", 0.0)? != 0.0,
            depth,
            shade: num("SHADE", 0.0)?,
            decay_chem: num("DECAY_CHEM", 0.0)?,
            decay_chem_seed: num("DECAY_CHEM_SEED", 1.0)? as u64,
            mutate,
            maxl: num("MAXLEN", 2000.0)? as usize,
            seed: num("RANDSEED", 1.0)? as u64,
            nsteps: num("NSTEPS", 1_000_000.0)? as u64,
            popdy_every: num("POPDYEVERY", 100.0)? as u64,
            snap_every: num("SNAPEVERY", 500.0)? as u64,
            copy_every: num("COPYEVERY", 1000.0)? as u64,
            ckpt_every: num("CKPTEVERY", 50_000.0)? as u64,
            keep_ckpt: num("KEEPCKPT", 3.0)? as usize,
            live_every: num("LIVEEVERY", 0.0)? as u64,
            bindlen: num("BINDLEN", 0.0)? as usize,
            bind_t: num("BIND_T", 0.0)?,
            compat_partner: num("COMPAT_PARTNER", 0.0)? != 0.0,
            decay_independent: num("DECAY_INDEPENDENT", 0.0)? != 0.0,
            encounter_turn: num("ENCOUNTER_TURN", 0.0)? != 0.0,
            displace: num("DISPLACE", 0.0)? != 0.0,
            displace_k: num("DISPLACE_K", 10.0)?,
            displace_rate: num("DISPLACE_RATE", 1.0)?,
            displace_restore: num("DISPLACE_RESTORE", 0.0)? != 0.0,
            twosided: num("TWOSIDED", 0.0)? != 0.0,
            strandcopy: num("STRANDCOPY", 0.0)? != 0.0,
            offspring_log: num("OFFSPRINGLOG", 0.0)? != 0.0,
            species_forget: num("SPECIESFORGET", 0.0)? != 0.0,
            disturb_every: num("DISTURB_EVERY", 0.0)? as u64,
            disturb_until: num("DISTURB_UNTIL", 0.0)? as u64,
            disturb_size: num("DISTURB_SIZE", 7.0)? as usize,
            disturb_mode: num("DISTURB_MODE", 1.0)? as u32,
            disturb_seed: num("DISTURB_SEED", 1.0)? as u64,
            gziplogs: num("GZIPLOGS", 0.0)? != 0.0,
            ckpt_keep_every: num("CKPTKEEPEVERY", 0.0)? as u64,
            translate: num("TRANSLATE", 0.0)? != 0.0,
            translate_seed: num("TRANSLATE_SEED", 1.0)? as u64,
            translate_stop: num("TRANSLATE_STOP", 0.0)? != 0.0,
            agents,
            raw,
        })
    }
}
