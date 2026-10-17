//! Alphabet, complements, the similarity table and mutation neighbours (spec: Symbols).

use crate::score_table::SCORE_TEXT;

/// The 33 symbols in their fixed circular order.
pub const KEY: &[u8; 33] = b"ABC$DEF%GH^IJK?LMN}OPQ>RST=UVWXYZ";
pub const NSYM: usize = 33;

pub struct Table {
    /// score[a][b] for symbol indexes a, b; score[NSYM][b] is the gap score for b.
    pub score: [[f32; NSYM]; NSYM + 1],
    /// byte -> symbol index, or -1
    pub index: [i8; 256],
}

impl Table {
    pub fn new() -> Self {
        let mut score = [[0f32; NSYM]; NSYM + 1];
        for (i, row) in SCORE_TEXT.iter().enumerate() {
            for (j, txt) in row.iter().enumerate() {
                score[i][j] = txt.parse::<f32>().expect("score table");
            }
        }
        let mut index = [-1i8; 256];
        for (i, &c) in KEY.iter().enumerate() {
            index[c as usize] = i as i8;
        }
        Table { score, index }
    }

    #[inline]
    pub fn idx(&self, c: u8) -> usize {
        let i = self.index[c as usize];
        debug_assert!(i >= 0, "unknown symbol {:?}", c as char);
        i as usize
    }

    /// Score for aligning symbol a with b (both non-zero).
    #[inline]
    pub fn pair(&self, a: u8, b: u8) -> f32 {
        self.score[self.idx(a)][self.idx(b)]
    }

    /// Score for a gap against symbol a.
    #[inline]
    pub fn gap(&self, a: u8) -> f32 {
        self.score[NSYM][self.idx(a)]
    }
}

impl Default for Table {
    fn default() -> Self {
        Self::new()
    }
}

#[inline]
pub fn is_letter(c: u8) -> bool {
    c.is_ascii_uppercase()
}

/// Complement: letters shift 13 places (A<->N), everything else is its own complement.
#[inline]
pub fn complement(c: u8) -> u8 {
    if is_letter(c) {
        (c - b'A' + 13) % 26 + b'A'
    } else {
        c
    }
}

/// The two neighbours of a symbol in the circular order, lower index first
/// (the order the old `sym_from_adj` scanned them in).
pub fn neighbours(c: u8, t: &Table) -> [u8; 2] {
    let i = t.idx(c);
    let a = KEY[(i + NSYM - 1) % NSYM];
    let b = KEY[(i + 1) % NSYM];
    let (ia, ib) = (t.idx(a), t.idx(b));
    if ia < ib { [a, b] } else { [b, a] }
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn complement_letters() {
        assert_eq!(complement(b'A'), b'N');
        assert_eq!(complement(b'N'), b'A');
        assert_eq!(complement(b'Z'), b'M');
        assert_eq!(complement(b'$'), b'$');
    }
    #[test]
    fn table_shape() {
        let t = Table::new();
        assert_eq!(t.pair(b'A', b'A'), 1.0);
        assert_eq!(t.pair(b'$', b'$'), 0.5);
        assert_eq!(t.pair(b'A', b'B'), -0.121);
        assert_eq!(t.gap(b'A'), -1.333);
    }
    #[test]
    fn neighbours_wrap() {
        let t = Table::new();
        assert_eq!(neighbours(b'A', &t), [b'B', b'Z']);
        assert_eq!(neighbours(b'C', &t), [b'B', b'$']);
    }
}
