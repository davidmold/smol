//! Codon table helper for TRANSLATE designs.
//!
//!   smol-codon encode STRING [SEED]   -> a gene: each symbol written as a letter pair
//!   smol-codon decode GENE [SEED]     -> what translating GENE writes
//!   smol-codon table [SEED]           -> each symbol and its codons
use smol::symbols::{codon_table, is_letter, KEY};

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let seed: u64 = a.get(3).or(if a.get(1).map(|s| s == "table").unwrap_or(false) { a.get(2) } else { None })
        .and_then(|v| v.parse().ok()).unwrap_or(1);
    let t = codon_table(seed);
    let codons = |c: u8| -> Vec<String> {
        (0..676).filter(|&i| t[i] == c).map(|i| format!("{}{}", (b'A' + (i / 26) as u8) as char, (b'A' + (i % 26) as u8) as char)).collect()
    };
    match a.get(1).map(|s| s.as_str()) {
        Some("encode") => {
            let s = a[2].as_bytes();
            let out: String = s.iter().map(|&c| codons(c)[0].clone()).collect();
            println!("{out}");
        }
        Some("decode") => {
            let g = a[2].as_bytes();
            let mut i = 0;
            let mut out = String::new();
            while i < g.len() {
                if i + 1 < g.len() && is_letter(g[i]) && is_letter(g[i + 1]) {
                    out.push(t[(g[i] - b'A') as usize * 26 + (g[i + 1] - b'A') as usize] as char);
                    i += 2;
                } else {
                    out.push(g[i] as char);
                    i += 1;
                }
            }
            println!("{out}");
        }
        _ => {
            for &c in KEY.iter() {
                println!("{} {}", c as char, codons(c).join(" "));
            }
        }
    }
}
