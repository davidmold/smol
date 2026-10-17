//! Output files (spec: Outputs and reproducibility).

use crate::world::{Bond, Observer, World};
use flate2::write::{GzEncoder, ZlibEncoder};
use flate2::Compression;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

pub struct Output {
    popdy: BufWriter<File>,
    species: BufWriter<File>,
    copies: BufWriter<File>,
    snap: BufWriter<File>,
    copylog: HashMap<(u32, u32, u32), u32>,
}

fn open(dir: &Path, name: &str) -> std::io::Result<BufWriter<File>> {
    Ok(BufWriter::new(OpenOptions::new().create(true).append(true).open(dir.join(name))?))
}

impl Output {
    pub fn new(dir: &Path) -> std::io::Result<Output> {
        Ok(Output {
            popdy: open(dir, "popdy001.dat")?,
            species: open(dir, "species.tsv")?,
            copies: open(dir, "copies.tsv")?,
            snap: open(dir, "grid.snap")?,
            copylog: HashMap::new(),
        })
    }

    pub fn popdy(&mut self, w: &World) -> std::io::Result<()> {
        let mut counts: HashMap<u32, u32> = HashMap::new();
        for m in w.mols.iter().filter(|m| m.alive) {
            *counts.entry(m.species).or_default() += 1;
        }
        let mut v: Vec<_> = counts.into_iter().collect();
        v.sort_unstable();
        for (sp, n) in v {
            writeln!(self.popdy, "{},{},{}", w.step, sp, n)?;
        }
        self.popdy.flush()
    }

    pub fn snapshot(&mut self, w: &World) -> std::io::Result<()> {
        let mut raw = Vec::with_capacity(w.grid.len() * 4);
        for &g in &w.grid {
            let v: u32 = if g == 0 {
                0
            } else {
                let m = &w.mols[(g - 1) as usize];
                m.species | if m.bond == Bond::Free { 0 } else { 0x8000_0000 }
            };
            raw.extend_from_slice(&v.to_le_bytes());
        }
        let mut z = ZlibEncoder::new(Vec::new(), Compression::default());
        z.write_all(&raw)?;
        let z = z.finish()?;
        for h in [w.step as u32, w.gx as u32, w.gy as u32, z.len() as u32] {
            self.snap.write_all(&h.to_le_bytes())?;
        }
        self.snap.write_all(&z)?;
        self.snap.flush()
    }

    pub fn flush_copies(&mut self, step: u64) -> std::io::Result<()> {
        let mut v: Vec<_> = self.copylog.drain().collect();
        v.sort_unstable();
        for ((a, p, c), n) in v {
            writeln!(self.copies, "{step}\t{a}\t{p}\t{c}\t{n}\tA")?;
        }
        self.copies.flush()?;
        self.species.flush()
    }
}

impl Observer for Output {
    fn new_species(&mut self, id: u32, born: u64, pa: i64, pp: i64, seq: &[u8]) {
        let _ = writeln!(self.species, "{}\t{}\t{}\t{}\t{}", id, born, pa, pp, String::from_utf8_lossy(seq));
    }
    fn product(&mut self, sp_act: u32, sp_pas: u32, sp_prod: u32) {
        *self.copylog.entry((sp_act, sp_pas, sp_prod)).or_default() += 1;
    }
}

/// Checkpoint: gzipped text. Species ids are restored from species.tsv on resume.
pub fn write_checkpoint(dir: &Path, w: &World, config_text: &str) -> std::io::Result<String> {
    let name = format!("checkpoint_{}.txt.gz", w.step);
    let f = File::create(dir.join(&name))?;
    let mut z = GzEncoder::new(BufWriter::new(f), Compression::default());
    writeln!(z, "SMOLCKPT 1")?;
    writeln!(z, "STEP {}", w.step)?;
    writeln!(z, "RNG {} {} {} {}", w.rng.s[0], w.rng.s[1], w.rng.s[2], w.rng.s[3])?;
    writeln!(z, "NEXTSPECIES {}", w.next_species)?;
    writeln!(z, "CONFIG {}", config_text.lines().count())?;
    for l in config_text.lines() {
        writeln!(z, "{l}")?;
    }
    writeln!(z, "SLOTS {}", w.mols.len())?;
    write!(z, "FREE {}", w.free_slots.len())?;
    for f in &w.free_slots {
        write!(z, " {f}")?;
    }
    writeln!(z)?;
    let live: Vec<usize> = (0..w.mols.len()).filter(|&i| w.mols[i].alive).collect();
    writeln!(z, "MOLS {}", live.len())?;
    for &i in &live {
        let m = &w.mols[i];
        // the whole buffer up to its last non-zero byte; zeros inside written as '.'
        let last = m.strand.buf.iter().rposition(|&c| c != 0).map(|p| p + 1).unwrap_or(0);
        let s: String = m.strand.buf[..last].iter().map(|&c| if c == 0 { '.' } else { c as char }).collect();
        let bond = match &m.bond {
            Bond::Free => "F".to_string(),
            Bond::Passive { partner } => format!("P {partner}"),
            Bond::Active { partner, ptr, sp_act, sp_pas, p0, copied } => format!(
                "A {partner} {sp_act} {sp_pas} {} {} {} {} {} {} {} {} {} {} {} {} {p0:e} {copied}",
                ptr.pos[0][0], ptr.pos[0][1], ptr.pos[1][0], ptr.pos[1][1],
                ptr.pos[2][0], ptr.pos[2][1], ptr.pos[3][0], ptr.pos[3][1],
                ptr.side[0], ptr.side[1], ptr.side[2], ptr.side[3]
            ),
        };
        writeln!(z, "{i} {} {} {} {s} {bond}", m.x, m.y, m.species)?;
    }
    z.finish()?.flush()?;
    Ok(name)
}
