//! Stringmol, reimplemented from the written specification.

pub mod align;
pub mod config;
pub mod output;
pub mod reaction;
pub mod rng;
mod score_table;
pub mod symbols;
pub mod world;

use align::{bind_prob, Align, Aligner};
use symbols::{complement, Table};

/// Align the complement of `a` against `b`, as when `a` tries to bind `b`.
pub fn bind_align(a: &[u8], b: &[u8], al: &mut Aligner, tab: &Table, comp: &mut Vec<u8>) -> Align {
    comp.clear();
    comp.extend(a.iter().map(|&c| complement(c)));
    al.align(comp, b, tab)
}

/// Binding probability of `a` initiating with `b`.
pub fn bind_probability(a: &[u8], b: &[u8], bindlen: usize, al: &mut Aligner, tab: &Table) -> (f32, Align) {
    let mut comp = Vec::with_capacity(a.len());
    let sw = bind_align(a, b, al, tab, &mut comp);
    (bind_prob(&sw, bindlen), sw)
}
