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
            let mut out = Output::new(dir).unwrap_or_else(|e| die(&e.to_string()));
            let w = World::new(&cfg, &mut out).unwrap_or_else(|e| die(&e));
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
            // species ids so far, from this directory's species.tsv
            if let Ok(sp) = std::fs::read_to_string(dir.join("species.tsv")) {
                for l in sp.lines() {
                    let p: Vec<&str> = l.split('\t').collect();
                    if p.len() >= 5 {
                        if let Ok(id) = p[0].parse::<u32>() {
                            if id < w.next_species {
                                w.known_species(id, p[4].as_bytes());
                            }
                        }
                    }
                }
            }
            let out = Output::new(dir).unwrap_or_else(|e| die(&e.to_string()));
            (cfg, cfg_text, w, out)
        }
        _ => die("usage: smol run CONF | smol resume CKPT [--reseed S] [--nsteps N]"),
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
        if t % cfg.copy_every == 0 {
            out.flush_copies(t).unwrap();
        }
        if t % cfg.ckpt_every == 0 {
            out.flush_copies(t).unwrap();
            ckpts.push(write_checkpoint(dir, &world, &cfg_text).unwrap());
            while ckpts.len() > cfg.keep_ckpt {
                let _ = std::fs::remove_file(dir.join(ckpts.remove(0)));
            }
        }
        if t % 1000 == 0 {
            let n = world.count();
            let bound = world.mols.iter().filter(|m| m.alive && m.bond != Bond::Free).count();
            println!("Step {t} done, number of molecules is {n}, nbound = {bound}, {:.0} steps/s", t as f64 / start.elapsed().as_secs_f64().max(1e-9));
        }
        if world.count() == 0 {
            out.popdy(&world).unwrap();
            println!("EXTINCT at step {t}");
            break;
        }
    }
    out.flush_copies(world.step).unwrap();
    write_checkpoint(dir, &world, &cfg_text).unwrap();
}
