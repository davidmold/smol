//! Where does a `$` search land? For each input line "LABEL TARGET" (LABEL = the letters after
//! `$` in the searching molecule, TARGET = the strand searched), prints the outcome exactly as the
//! `$` instruction computes it (reaction.rs): the complemented label is aligned to the target
//! (Smith-Waterman), and with probability `prob` the flow pointer moves to `pos` (the end of the
//! alignment in the target); otherwise it stays on the searching strand.
//!
//!   smol-search < lines      ->  "pos prob score aligned_length" per line
use smol::align::{search_prob, Aligner};
use smol::symbols::{complement, Table};
use std::io::BufRead;

fn main() {
    let tab = Table::new();
    let mut al = Aligner::new();
    for line in std::io::stdin().lock().lines() {
        let line = line.unwrap();
        let mut it = line.split_whitespace();
        let (Some(label), Some(target)) = (it.next(), it.next()) else { continue };
        let lab: Vec<u8> = label.bytes().map(complement).collect();
        let a = al.align(&lab, target.as_bytes(), &tab);
        println!("{} {:.4} {:.3} {}", a.e2, search_prob(&a), a.score, a.len());
    }
}
