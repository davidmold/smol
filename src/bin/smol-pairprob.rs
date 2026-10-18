//! Binding probability for each line "SEQ_A SEQ_B" on stdin (A initiating with B), plus the
//! alignment: prob score s1 e1 s2 e2. For checking against the old code.
use smol::align::Aligner;
use smol::symbols::Table;
use std::io::{BufRead, Write};

fn main() {
    let bindlen: usize = std::env::args().nth(1).and_then(|s| s.parse().ok()).unwrap_or(0);
    let bind_t: f64 = std::env::args().nth(2).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let tab = Table::new();
    let mut al = Aligner::new();
    let stdin = std::io::stdin();
    let out = std::io::stdout();
    let mut out = out.lock();
    for line in stdin.lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let (Some(a), Some(b)) = (it.next(), it.next()) else { continue };
        let (p, sw) = smol::bind_probability_rule(a.as_bytes(), b.as_bytes(), bindlen, bind_t, &mut al, &tab);
        writeln!(out, "{:.4} {:.4} {} {} {} {}", p, sw.score, sw.s1, sw.e1, sw.s2, sw.e2).unwrap();
    }
}
