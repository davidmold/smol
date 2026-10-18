//! Reaction probe: react every listed molecule with every other (and itself), each taking a
//! turn as initiator, mutation off, and report who copied whom and what was produced.
//! Output format matches the old `sm-react`, so the two can be diffed.
//!
//!   smol-react PAIRS [MAXSTEPS] [--compat] [--skip N]
//!
//! PAIRS lines: "R name seq" or "P name seq". --compat uses the old generator (MT19937,
//! default seed 4357) and the old deletion quirk, for draw-for-draw comparison; --skip N
//! discards N draws first, to line up with draws the old code spent while loading.

use smol::align::Aligner;
use smol::reaction::{clamp_ptrs, exec_step, ChemParams, Ctx, Outcome, Pointers, Strand};
use smol::rng::{Mt19937, Rng, Xoshiro};
use smol::symbols::Table;
use std::io::Write;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        eprintln!("usage: smol-react PAIRS [MAXSTEPS] [--compat] [--skip N]");
        std::process::exit(1);
    }
    let compat = args.iter().any(|a| a == "--compat");
    let skip: usize = args
        .iter()
        .position(|a| a == "--skip")
        .map(|i| args[i + 1].parse().unwrap())
        .unwrap_or(0);
    let maxsteps: usize = args.get(2).and_then(|s| s.parse().ok()).unwrap_or(20000);
    let mutate: f32 = args
        .iter()
        .position(|a| a == "--mutate")
        .map(|i| args[i + 1].parse().unwrap())
        .unwrap_or(0.0);
    let bind_t: f64 = args.iter().position(|a| a == "--bind-t").map(|i| args[i + 1].parse().unwrap()).unwrap_or(0.0);
    let twosided = args.iter().any(|a| a == "--twosided");
    let trace = args.iter().any(|a| a == "--trace");
    let text = std::fs::read_to_string(&args[1]).expect("pairs file");
    let mols: Vec<(String, Vec<u8>)> = text
        .split_whitespace()
        .collect::<Vec<_>>()
        .chunks(3)
        .filter(|c| c.len() == 3)
        .map(|c| (c[1].to_string(), c[2].as_bytes().to_vec()))
        .collect();

    let mut rng: Box<dyn Rng> = if compat { Box::new(Mt19937::new(4357)) } else { Box::new(Xoshiro::new(1)) };
    for _ in 0..skip {
        rng.next_f64();
    }
    let tab = Table::new();
    let mut al = Aligner::new();
    let params = ChemParams { maxl: 2000, mutate, compat_deletion_skip: compat };
    let mut ctx = Ctx { tab: &tab, al: &mut al, params, label: Vec::new(), other: Vec::new() };
    let out = std::io::stdout();
    let mut out = out.lock();

    let ident = |s: &[u8]| -> String {
        mols.iter()
            .find(|m| m.1 == s)
            .map(|m| m.0.clone())
            .unwrap_or_else(|| format!("new({})", s.len()))
    };

    for x in &mols {
        for y in &mols {
            let (p, sw) = smol::bind_probability_rule(&x.1, &y.1, 0, bind_t, ctx.al, &tab);
            // old set_exec: the molecule whose alignment starts later is active
            let (act_is_x, act_idx, pas_idx) = if sw.s1 >= sw.s2 { (true, sw.s1, sw.s2) } else { (false, sw.s2, sw.s1) };
            let (an, pn) = if act_is_x { (&x.0, &y.0) } else { (&y.0, &x.0) };
            let mut own = Strand::new(if act_is_x { &x.1 } else { &y.1 }, params.maxl);
            let mut partner = Strand::new(if act_is_x { &y.1 } else { &x.1 }, params.maxl);
            let mut ptr = Pointers::at_bind(act_idx, pas_idx);
            let mut pptr = Pointers::at_bind(pas_idx, act_idx);
            let mut products = Vec::new();
            let mut pprods: Vec<Vec<u8>> = Vec::new();
            let mut steps = 0;
            let mut end = "ran out of steps";
            let mut end_by_pas = false;
            while steps < maxsteps {
                let mut o = exec_step(&mut ptr, &mut own, &mut partner, &mut ctx, rng.as_mut(), &mut products);
                if twosided && o == Outcome::Continue {
                    clamp_ptrs(&mut pptr, partner.len(), own.len());
                    o = match exec_step(&mut pptr, &mut partner, &mut own, &mut ctx, rng.as_mut(), &mut pprods) {
                        Outcome::Ended => { end_by_pas = true; Outcome::Ended }
                        Outcome::ActiveDestroyed => Outcome::PassiveDestroyed,
                        Outcome::PassiveDestroyed => Outcome::ActiveDestroyed,
                        c => c,
                    };
                    clamp_ptrs(&mut ptr, own.len(), partner.len());
                }
                steps += 1;
                if trace && (steps <= 400 && steps % 20 == 0 || steps % 1000 == 0) {
                    eprintln!("  t={steps} act len {} I{:?} R{:?} W{:?} | pas len {} I{:?} R{:?} W{:?}", own.len(),
                        (ptr.side[0], ptr.pos[0][ptr.side[0]]), (ptr.side[1], ptr.pos[1][ptr.side[1]]), (ptr.side[2], ptr.pos[2][ptr.side[2]]),
                        partner.len(), (pptr.side[0], pptr.pos[0][pptr.side[0]]), (pptr.side[1], pptr.pos[1][pptr.side[1]]), (pptr.side[2], pptr.pos[2][pptr.side[2]]));
                }
                match o {
                    Outcome::Continue => {}
                    Outcome::Ended => {
                        end = "unbound normally";
                        break;
                    }
                    Outcome::ActiveDestroyed => {
                        end = "active destroyed";
                        break;
                    }
                    Outcome::PassiveDestroyed => {
                        end = "passive destroyed";
                        break;
                    }
                }
            }
            writeln!(
                out,
                "{:<8} initiates with {:<8} P(bind)={:.3}  copier={:<8} copied={:<8} steps={:<5} {}",
                x.0, y.0, p, an, pn, steps, end
            )
            .unwrap();
            if end_by_pas {
                writeln!(out, "      (ended by the passive)").unwrap();
            }
            for s in &products {
                writeln!(out, "      -> product {:<10} {}", ident(s), String::from_utf8_lossy(s)).unwrap();
            }
            for s in &pprods {
                writeln!(out, "      -> passive made {:<10} {}", ident(s), String::from_utf8_lossy(s)).unwrap();
            }
            if products.is_empty() && pprods.is_empty() {
                writeln!(out, "      -> no product").unwrap();
            }
        }
    }
}
