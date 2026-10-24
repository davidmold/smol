//! Output files (spec: Outputs and reproducibility).

use crate::world::{Bond, Observer, World};
use flate2::write::{GzEncoder, ZlibEncoder};
use flate2::Compression;
use std::collections::HashMap;
use std::fs::{File, OpenOptions};
use std::io::{BufWriter, Write};
use std::path::Path;

/// A log file: plain text, or (GZIPLOGS) gzip written as a new member at every checkpoint, so a
/// run that is killed loses at most the data since its last checkpoint.
pub trait Sink: Write {
    fn rotate(&mut self) -> std::io::Result<()> {
        self.flush()
    }
}
impl Sink for BufWriter<File> {}

pub struct GzSink {
    enc: Option<GzEncoder<BufWriter<File>>>,
    path: std::path::PathBuf,
}
impl GzSink {
    fn open(path: std::path::PathBuf) -> std::io::Result<GzSink> {
        let f = BufWriter::new(OpenOptions::new().create(true).append(true).open(&path)?);
        Ok(GzSink { enc: Some(GzEncoder::new(f, Compression::default())), path })
    }
}
impl Write for GzSink {
    fn write(&mut self, b: &[u8]) -> std::io::Result<usize> {
        self.enc.as_mut().unwrap().write(b)
    }
    fn flush(&mut self) -> std::io::Result<()> {
        Ok(()) // a gzip member is completed only at rotate() or drop
    }
}
impl Sink for GzSink {
    fn rotate(&mut self) -> std::io::Result<()> {
        if let Some(e) = self.enc.take() {
            e.finish()?.flush()?;
        }
        let f = BufWriter::new(OpenOptions::new().create(true).append(true).open(&self.path)?);
        self.enc = Some(GzEncoder::new(f, Compression::default()));
        Ok(())
    }
}
impl Drop for GzSink {
    fn drop(&mut self) {
        if let Some(e) = self.enc.take() {
            if let Ok(mut f) = e.finish() {
                let _ = f.flush();
            }
        }
    }
}

pub struct Output {
    popdy: Box<dyn Sink>,
    species: Box<dyn Sink>,
    copies: Box<dyn Sink>,
    snap: BufWriter<File>,
    copylog: HashMap<(u32, u32, u32), (u32, u32)>,
    gz: bool,
    /// OFFSPRINGLOG: write the placed column and lives.tsv
    pub offspring: Option<Box<dyn Sink>>,
}

fn open(dir: &Path, name: &str) -> std::io::Result<BufWriter<File>> {
    Ok(BufWriter::new(OpenOptions::new().create(true).append(true).open(dir.join(name))?))
}

fn sink(dir: &Path, name: &str, gz: bool) -> std::io::Result<Box<dyn Sink>> {
    Ok(if gz { Box::new(GzSink::open(dir.join(format!("{name}.gz")))?) } else { Box::new(open(dir, name)?) })
}

impl Output {
    /// gz (GZIPLOGS): popdy001.dat, species.tsv, copies.tsv and lives.tsv are written gzipped.
    pub fn new(dir: &Path, gz: bool) -> std::io::Result<Output> {
        Ok(Output {
            popdy: sink(dir, "popdy001.dat", gz)?,
            species: sink(dir, "species.tsv", gz)?,
            copies: sink(dir, "copies.tsv", gz)?,
            snap: open(dir, "grid.snap")?,
            copylog: HashMap::new(),
            gz,
            offspring: None,
        })
    }

    /// Close the current gzip members (GZIPLOGS); called after each checkpoint.
    pub fn rotate_logs(&mut self) -> std::io::Result<()> {
        self.popdy.rotate()?;
        self.species.rotate()?;
        self.copies.rotate()?;
        if let Some(f) = &mut self.offspring {
            f.rotate()?;
        }
        Ok(())
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

    /// Latest state for a live viewer, replaced atomically: live.grid (one snapshot frame, same
    /// layout as a grid.snap frame) and live.json (step, counts, the 40 most abundant species).
    pub fn live(&mut self, dir: &Path, w: &World) -> std::io::Result<()> {
        // DEPTH > 1: slices are drawn as panels stacked top to bottom, one empty row between
        let (pw, ph) = if w.gz > 1 { (w.sx, w.gz * w.gy + w.gz - 1) } else { (w.gx, w.gy) };
        let mut raw = Vec::with_capacity(pw * ph * 4);
        for r in 0..ph {
            let (z, y) = (r / (w.gy + 1), r % (w.gy + 1));
            for c in 0..pw {
                let g = if w.gz > 1 {
                    if y == w.gy { 0 } else { w.grid[y * w.gx + z * w.sx + c] }
                } else {
                    w.grid[r * w.gx + c]
                };
                let v: u32 = if g == 0 {
                    0
                } else {
                    let m = &w.mols[(g - 1) as usize];
                    m.species | if m.bond == Bond::Free { 0 } else { 0x8000_0000 }
                };
                raw.extend_from_slice(&v.to_le_bytes());
            }
        }
        let mut z = ZlibEncoder::new(Vec::new(), Compression::fast());
        z.write_all(&raw)?;
        let z = z.finish()?;
        let mut buf = Vec::with_capacity(z.len() + 16);
        for h in [w.step as u32, pw as u32, ph as u32, z.len() as u32] {
            buf.extend_from_slice(&h.to_le_bytes());
        }
        buf.extend_from_slice(&z);
        std::fs::write(dir.join("live.grid.tmp"), &buf)?;
        std::fs::rename(dir.join("live.grid.tmp"), dir.join("live.grid"))?;

        let mut counts: HashMap<u32, (u32, usize)> = HashMap::new();
        let (mut n, mut bound) = (0u32, 0u32);
        for (i, m) in w.mols.iter().enumerate() {
            if m.alive {
                n += 1;
                if m.bond != Bond::Free {
                    bound += 1;
                }
                counts.entry(m.species).or_insert((0, i)).0 += 1;
            }
        }
        let mut v: Vec<_> = counts.into_iter().collect();
        v.sort_unstable_by(|a, b| b.1 .0.cmp(&a.1 .0).then(a.0.cmp(&b.0)));
        let nspecies = v.len();
        let top: Vec<String> = v
            .iter()
            .take(40)
            .map(|(sp, (c, i))| {
                let seq = String::from_utf8_lossy(w.mols[*i].strand.seq()).into_owned();
                format!("{{\"id\":{sp},\"n\":{c},\"seq\":\"{seq}\"}}")
            })
            .collect();
        let json = format!(
            "{{\"step\":{},\"gx\":{},\"gy\":{},\"n\":{n},\"bound\":{bound},\"species\":{nspecies},\"top\":[{}]}}",
            w.step, pw, ph, top.join(",")
        );
        std::fs::write(dir.join("live.json.tmp"), json)?;
        std::fs::rename(dir.join("live.json.tmp"), dir.join("live.json"))
    }

    /// Turn on OFFSPRINGLOG output (lives.tsv in dir).
    pub fn enable_offspring(&mut self, dir: &Path) -> std::io::Result<()> {
        self.offspring = Some(sink(dir, "lives.tsv", self.gz)?);
        Ok(())
    }

    /// OFFSPRINGLOG: one row per species with deaths since the last call:
    /// step, species, died, mean lifetime, died having made >= 1 placed product, placed products
    /// made, died having been copied exactly >= 1 time (placed), placed exact copies.
    pub fn flush_lives(&mut self, w: &mut World) -> std::io::Result<()> {
        if let Some(f) = &mut self.offspring {
            let mut v: Vec<_> = w.lives.drain().collect();
            v.sort_unstable_by_key(|(k, _)| *k);
            for (sp, e) in v {
                writeln!(f, "{}\t{}\t{}\t{:.1}\t{}\t{}\t{}\t{}", w.step, sp, e.died, e.life_sum as f64 / e.died.max(1) as f64, e.made_any, e.made_sum, e.copied_any, e.copied_sum)?;
            }
            f.flush()?;
        }
        Ok(())
    }

    pub fn flush_copies(&mut self, step: u64) -> std::io::Result<()> {
        let mut v: Vec<_> = self.copylog.drain().collect();
        v.sort_unstable();
        let placed_col = self.offspring.is_some();
        for ((a, p, c), (n, k)) in v {
            if placed_col {
                writeln!(self.copies, "{step}\t{a}\t{p}\t{c}\t{n}\tA\t{k}")?;
            } else {
                writeln!(self.copies, "{step}\t{a}\t{p}\t{c}\t{n}\tA")?;
            }
        }
        self.copies.flush()?;
        self.species.flush()
    }
}

impl Observer for Output {
    fn new_species(&mut self, id: u32, born: u64, pa: i64, pp: i64, seq: &[u8]) {
        let _ = writeln!(self.species, "{}\t{}\t{}\t{}\t{}", id, born, pa, pp, String::from_utf8_lossy(seq));
    }
    fn product(&mut self, sp_act: u32, sp_pas: u32, sp_prod: u32, placed: bool) {
        let e = self.copylog.entry((sp_act, sp_pas, sp_prod)).or_default();
        e.0 += 1;
        e.1 += placed as u32;
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
            Bond::Active { partner, ptr, sp_act, sp_pas, p0, copied, pptr, nas, pnas } => format!(
                "A {partner} {sp_act} {sp_pas} {} {} {} {} {} {} {} {} {} {} {} {} {p0:e} {copied}",
                ptr.pos[0][0], ptr.pos[0][1], ptr.pos[1][0], ptr.pos[1][1],
                ptr.pos[2][0], ptr.pos[2][1], ptr.pos[3][0], ptr.pos[3][1],
                ptr.side[0], ptr.side[1], ptr.side[2], ptr.side[3]
            ) + &match pptr {
                None => String::new(),
                Some(q) => format!(
                    " {} {} {} {} {} {} {} {} {} {} {} {}",
                    q.pos[0][0], q.pos[0][1], q.pos[1][0], q.pos[1][1],
                    q.pos[2][0], q.pos[2][1], q.pos[3][0], q.pos[3][1],
                    q.side[0], q.side[1], q.side[2], q.side[3]
                ),
            } + &if w.strandcopy {
                let t = |v: &Vec<u8>| if v.is_empty() { "-".to_string() } else { String::from_utf8_lossy(v).into_owned() };
                format!(" n={} pn={}", t(nas), t(pnas))
            } else {
                String::new()
            },
        };
        writeln!(z, "{i} {} {} {} {s} {bond}", m.x, m.y, m.species)?;
    }
    z.finish()?.flush()?;
    Ok(name)
}
