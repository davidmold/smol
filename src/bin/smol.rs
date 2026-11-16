//! smol run CONF            start a run in the current directory
//! smol resume CKPT [--reseed S] [--nsteps N]
//!                          continue from a checkpoint in the current directory; --reseed
//!                          replaces the random state with a fresh seed (replay experiments)

use smol::config::Config;
use smol::output::{write_checkpoint, Output};
use smol::rng::Xoshiro;
use smol::world::{Bond, World};
use std::io::Read;
use std::path::Path;

fn die(msg: &str) -> ! {
    eprintln!("smol: {msg}");
    std::process::exit(1)
}

fn flag(args: &[String], name: &str) -> Option<u64> {
    args.iter().position(|a| a == name).and_then(|i| args.get(i + 1)).map(|v| v.parse().unwrap_or_else(|_| die(&format!("bad {name}"))))
}

fn main() {
    let args: Vec<String> = std::env::args().collect();
    let dir = Path::new(".");
    let (cfg, cfg_text, mut world, mut out) = match args.get(1).map(String::as_str) {
        Some("run") => {
            let path = args.get(2).unwrap_or_else(|| die("usage: smol run CONF"));
            let text = std::fs::read_to_string(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
            let cfg = Config::parse(&text).unwrap_or_else(|e| die(&e));
            for f in ["popdy001.dat", "species.tsv", "copies.tsv", "grid.snap"] {
                if dir.join(f).exists() {
                    die(&format!("{f} already exists here; run in an empty directory"));
                }
            }
            let mut out = Output::new(dir, cfg.gziplogs).unwrap_or_else(|e| die(&e.to_string()));
            let w = World::new(&cfg, &mut out).unwrap_or_else(|e| die(&e));
            if let Some(dw) = &w.dweights {
                let mut t = String::from("symbol\tweight\n");
                for &c in smol::symbols::KEY.iter() {
                    t += &format!("{}\t{:.4}\n", c as char, dw[c as usize]);
                }
                let _ = std::fs::write("decay_weights.tsv", t);
            }
            out.snapshot(&w).unwrap();
            (cfg, text, w, out)
        }
        Some("resume") => {
            let path = args.get(2).unwrap_or_else(|| die("usage: smol resume CKPT"));
            let f = std::fs::File::open(path).unwrap_or_else(|e| die(&format!("{path}: {e}")));
            let mut text = String::new();
            flate2::read::GzDecoder::new(f).read_to_string(&mut text).unwrap_or_else(|e| die(&e.to_string()));
            // the checkpoint carries its config
            let lines: Vec<&str> = text.lines().collect();
            let ci = lines.iter().position(|l| l.starts_with("CONFIG ")).unwrap_or_else(|| die("no CONFIG in checkpoint"));
            let n: usize = lines[ci][7..].trim().parse().unwrap();
            let cfg_text = lines[ci + 1..ci + 1 + n].join("\n");
            let mut cfg = Config::parse(&cfg_text).unwrap_or_else(|e| die(&e));
            if let Some(ns) = flag(&args, "--nsteps") {
                cfg.nsteps = ns;
            }
            let mut w = World::restore(&cfg, &text).unwrap_or_else(|e| die(&e));
            if let Some(s) = flag(&args, "--reseed") {
                w.rng = Xoshiro::new(s);
            }
            // species ids so far, from this directory's species.tsv and species.tsv.gz;
            // with SPECIESFORGET only the species now present (a molecule's strand can differ
            // from its species' sequence mid-reaction, so the sequences come from the file)
            {
                let live = if cfg.species_forget { Some(w.live_species()) } else { None };
                use std::io::BufRead;
                let mut readers: Vec<Box<dyn BufRead>> = Vec::new();
                if let Ok(f) = std::fs::File::open(dir.join("species.tsv")) {
                    readers.push(Box::new(std::io::BufReader::new(f)));
                }
                if let Ok(f) = std::fs::File::open(dir.join("species.tsv.gz")) {
                    readers.push(Box::new(std::io::BufReader::new(flate2::read::MultiGzDecoder::new(f))));
                }
                for r in readers {
                    for l in r.lines() {
                        let Ok(l) = l else { break }; // a truncated final gzip member ends the read
                        let p: Vec<&str> = l.split('\t').collect();
                        if p.len() >= 5 {
                            if let Ok(id) = p[0].parse::<u32>() {
                                if id < w.next_species && live.as_ref().map_or(true, |s| s.contains(&id)) {
                                    w.known_species(id, p[4].as_bytes());
                                }
                            }
                        }
                    }
                }
                w.forget_cutoff = w.next_species;
            }
            let out = Output::new(dir, cfg.gziplogs).unwrap_or_else(|e| die(&e.to_string()));
            (cfg, cfg_text, w, out)
        }
        _ => die("usage: smol run CONF | smol resume CKPT [--reseed S] [--nsteps N]"),
    };
    if cfg.offspring_log {
        out.enable_offspring(dir).unwrap_or_else(|e| die(&e.to_string()));
    }
    // TRACESEQ SEQ [TRACEMAX N]: log binds, reaction ends, products and deaths of the molecules
    // carrying SEQ now, and of their descendants, to trace.tsv (see docs/switches.md)
    let rawv = |k: &str| cfg.raw.iter().rev().find(|(key, _)| key == k).map(|(_, v)| v.trim().to_string());
    if let Some(ts) = rawv("TRACESEQ") {
        let max = rawv("TRACEMAX").and_then(|v| v.parse::<f64>().ok()).unwrap_or(200000.0) as usize;
        let at: Option<Vec<(usize, usize)>> = rawv("TRACEAT").map(|v| {
            v.split(';').filter_map(|c| {
                let mut it = c.split(',').map(|n| n.trim().parse::<usize>());
                match (it.next(), it.next()) { (Some(Ok(x)), Some(Ok(y))) => Some((x, y)), _ => None }
            }).collect()
        });
        world.tag_seq(ts.as_bytes(), max, at.as_deref());
    }
    let trace_out = |w: &mut World| {
        if let Some(t) = &mut w.trace {
            if !t.buf.is_empty() {
                use std::io::Write;
                let mut f = std::fs::OpenOptions::new().create(true).append(true).open(dir.join("trace.tsv")).unwrap();
                for l in t.buf.drain(..) {
                    let _ = writeln!(f, "{l}");
                }
            }
        }
    };

    let mut ckpts: Vec<String> = Vec::new();
    let start = std::time::Instant::now();
    while world.step < cfg.nsteps {
        world.step(&mut out);
        let t = world.step;
        if t % cfg.popdy_every == 0 {
            out.popdy(&world).unwrap();
        }
        if t % cfg.snap_every == 0 {
            out.snapshot(&world).unwrap();
        }
        if cfg.live_every > 0 && t % cfg.live_every == 0 {
            out.live(dir, &world).unwrap();
        }
        if t % cfg.copy_every == 0 {
            out.flush_copies(t).unwrap();
            out.flush_lives(&mut world).unwrap();
            trace_out(&mut world);
        }
        if t % cfg.ckpt_every == 0 {
            out.flush_copies(t).unwrap();
            ckpts.push(write_checkpoint(dir, &world, &cfg_text).unwrap());
            out.flush_lives(&mut world).unwrap();
            out.rotate_logs().unwrap();
            if cfg.species_forget {
                let n = world.forget_species();
                println!("SPECIESFORGET step {t}: forgot {n}, kept {}", world.species.len());
            }
            while ckpts.len() > cfg.keep_ckpt {
                let old = ckpts.remove(0);
                let step_of = |n: &str| n.trim_start_matches("checkpoint_").trim_end_matches(".txt.gz").parse::<u64>().unwrap_or(1);
                if cfg.ckpt_keep_every > 0 && step_of(&old) % cfg.ckpt_keep_every == 0 {
                    continue; // CKPTKEEPEVERY: a permanent checkpoint
                }
                let _ = std::fs::remove_file(dir.join(old));
            }
        }
        if t % 1000 == 0 {
            let n = world.count();
            let bound = world.mols.iter().filter(|m| m.alive && m.bond != Bond::Free).count();
            let st = std::mem::take(&mut world.stats);
            println!(
                "Step {t} done, number of molecules is {n}, nbound = {bound}, {:.0} steps/s | enc {} bind {} kin {} ended {} prod {} pairdecay {} displace {} instr {} damaged {} pprod {} pend {} dmg {} blocked {}",
                t as f64 / start.elapsed().as_secs_f64().max(1e-9),
                st.encounters, st.binds, st.kin_binds, st.ended, st.products, st.pair_decays, st.displacements, st.instr, st.disp_damaged, st.pas_products, st.pas_ended, st.damage_hits, st.blocked
            );
            if cfg.disturb_every > 0 && st.disturbed > 0 {
                println!("DISTURB to step {t}: removed {}", st.disturbed);
            }
        }
        if world.count() == 0 {
            out.popdy(&world).unwrap();
            println!("EXTINCT at step {t}");
            break;
        }
    }
    out.flush_copies(world.step).unwrap();
    out.flush_lives(&mut world).unwrap();
    trace_out(&mut world);
    write_checkpoint(dir, &world, &cfg_text).unwrap();
}
