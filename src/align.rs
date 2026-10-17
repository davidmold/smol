//! Local (Smith–Waterman) alignment, as the old `SmithWatermanV2`, and the probability
//! rules built on it. Arithmetic is f32 where the old code used float, and f64 where its
//! expressions were promoted to double, so results match it bit for bit.

use crate::symbols::Table;

#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Align {
    pub score: f32,
    /// DP indexes: the aligned region of string 1 is s1..e1 (0-based, end exclusive).
    pub s1: usize,
    pub e1: usize,
    pub s2: usize,
    pub e2: usize,
    pub matches: usize,
}

impl Align {
    /// Aligned length: the shorter of the two aligned spans.
    pub fn len(&self) -> usize {
        (self.e1 - self.s1).min(self.e2 - self.s2)
    }
}

const MATCH: u8 = 1;
const DEL: u8 = 2;
const INS: u8 = 3;

/// Reusable scratch space, so aligning allocates nothing in steady state.
#[derive(Default)]
pub struct Aligner {
    h: Vec<f32>,
    t: Vec<u8>,
}

impl Aligner {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn align(&mut self, s1: &[u8], s2: &[u8], tab: &Table) -> Align {
        let (l1, l2) = (s1.len(), s2.len());
        if l1 == 0 || l2 == 0 {
            // Spec quirk 4: an empty string matches nothing, explicitly.
            return Align::default();
        }
        let w = l2 + 1;
        let n = (l1 + 1) * w;
        self.h.clear();
        self.h.resize(n, 0.0);
        self.t.clear();
        self.t.resize(n, 0);
        let (h, t) = (&mut self.h, &mut self.t);
        let mut emax = 0f32;
        let (mut endi, mut endj) = (0usize, 0usize);
        for i in 1..=l1 {
            let a = s1[i - 1];
            let ga = tab.gap(a);
            let row = tab.idx(a);
            for j in 1..=l2 {
                let b = s2[j - 1];
                let k = i * w + j;
                let mut hv = 0f32;
                let mut tv = 0u8;
                let m = h[(i - 1) * w + j - 1] + tab.score[row][tab.idx(b)];
                if m > hv {
                    hv = m;
                    tv = MATCH;
                }
                let m = h[(i - 1) * w + j] + ga;
                if m > hv {
                    hv = m;
                    tv = DEL;
                }
                let m = h[k - 1] + tab.gap(b);
                if m > hv {
                    hv = m;
                    tv = INS;
                }
                h[k] = hv;
                t[k] = tv;
                if hv > emax {
                    emax = hv;
                    endi = i;
                    endj = j;
                }
            }
        }
        let (mut i, mut j, mut matches) = (endi, endj, 0usize);
        loop {
            match t[i * w + j] {
                MATCH => {
                    matches += 1;
                    i -= 1;
                    j -= 1;
                }
                DEL => i -= 1,
                INS => j -= 1,
                _ => break,
            }
        }
        Align { score: emax, s1: i, e1: endi, s2: j, e2: endj, matches }
    }
}

/// Binding probability (spec: Binding, step 3). `bindlen` 0 is the old rule.
pub fn bind_prob(al: &Align, bindlen: usize) -> f32 {
    let l = al.len();
    if l <= 2 {
        return 0.0;
    }
    let lb = l.max(bindlen) as f64 - 1.124;
    let s: f32 = if (al.score as f64) < lb { al.score } else { lb as f32 };
    (s as f64 / lb) as f32
}

/// Success probability of a `$` search. The old code computes the binding formula but
/// without the l <= 2 cut-off (its zeroing is overwritten), so short matches succeed
/// almost surely.
pub fn search_prob(al: &Align) -> f32 {
    let lb = al.len() as f64 - 1.124;
    let s: f32 = if (al.score as f64) < lb { al.score } else { lb as f32 };
    (s as f64 / lb) as f32
}

/// Success probability of a `?` label test of length `len` (old `align_event`).
pub fn label_prob(al: &Align, len: usize) -> f32 {
    let x = al.score / len as f32;
    (x as f64).powf(len as f64) as f32
}
