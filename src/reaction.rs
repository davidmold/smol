//! One reaction: pointer state and the opcode interpreter (spec: Execution, Copying and
//! mutation). Strings are fixed buffers of `maxl + 1` bytes, zero-filled, exactly as in the
//! old code, so overwrites past the end, gaps and truncation behave identically. Pointers are
//! indexes into those buffers, never addresses.

use crate::align::{label_prob, search_prob, Aligner};
use crate::rng::Rng;
use crate::symbols::{complement, is_letter, neighbours, Table, KEY, NSYM};

pub const I: usize = 0;
pub const R: usize = 1;
pub const W: usize = 2;
pub const F: usize = 3;

/// Side 1 is the executing (active) molecule's own string, side 0 its partner's.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Pointers {
    pub pos: [[usize; 2]; 4],
    pub side: [usize; 4],
}

impl Pointers {
    /// Spec: Binding, step 5.
    pub fn at_bind(own_start: usize, partner_start: usize) -> Self {
        Pointers { pos: [[partner_start, own_start]; 4], side: [1; 4] }
    }
    #[inline]
    fn get(&self, p: usize) -> usize {
        self.pos[p][self.side[p]]
    }
    #[inline]
    fn set(&mut self, p: usize, v: usize) {
        self.pos[p][self.side[p]] = v;
    }
    #[inline]
    fn bump(&mut self, p: usize, by: usize) {
        self.pos[p][self.side[p]] += by;
    }
}

/// A molecule's string buffer.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Strand {
    pub buf: Vec<u8>,
}

impl Strand {
    pub fn new(seq: &[u8], maxl: usize) -> Self {
        let mut buf = vec![0u8; maxl + 1];
        let n = seq.len().min(maxl);
        buf[..n].copy_from_slice(&seq[..n]);
        Strand { buf }
    }
    #[inline]
    pub fn at(&self, i: usize) -> u8 {
        if i < self.buf.len() { self.buf[i] } else { 0 }
    }
    #[inline]
    pub fn put(&mut self, i: usize, c: u8) {
        if i < self.buf.len() {
            self.buf[i] = c;
        }
    }
    /// C strlen: up to the first zero.
    pub fn len(&self) -> usize {
        self.buf.iter().position(|&c| c == 0).unwrap_or(self.buf.len())
    }
    pub fn seq(&self) -> &[u8] {
        &self.buf[..self.len()]
    }
    pub fn is_empty(&self) -> bool {
        self.buf[0] == 0
    }
}

#[derive(Clone, Copy, Debug)]
pub struct ChemParams {
    pub maxl: usize,
    /// Probability of an indel per copied symbol, and separately of a substitution.
    pub mutate: f32,
    /// Old behaviour for validation only: a deletion also skips the next instruction.
    pub compat_deletion_skip: bool,
    /// TRANSLATE: the codon table (letter pair -> symbol). Only used with STRANDCOPY, when the
    /// write pointer's side has been toggled away from the executor's own string.
    pub codon: Option<[u8; 676]>,
    /// TRANSLATE_STOP: in translate mode, a symbol under R that is not the start of a letter
    /// pair is a stop codon: the reaction ends and the nascent strand is released.
    pub translate_stop: bool,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Outcome {
    Continue,
    /// The reaction ended (`}`, end of code, or a copy past the length limit).
    Ended,
    /// A cleave left the active's string empty: it is destroyed, the passive freed.
    ActiveDestroyed,
    /// A cleave left the passive's string empty.
    PassiveDestroyed,
}

/// Scratch shared by all reactions.
pub struct Ctx<'a> {
    pub tab: &'a Table,
    pub al: &'a mut Aligner,
    pub params: ChemParams,
    pub label: Vec<u8>,
    pub other: Vec<u8>,
    /// STRANDCOPY: the executing molecule's nascent copy strand. `=` appends to it instead of
    /// writing at W, and `%` releases it as a product. None = the old in-place copying.
    pub nascent: Option<Vec<u8>>,
}

/// After the other side of a two-sided reaction changed a string, pull this side's pointers
/// that now point past an end back to it (the old rewind_bad_ptrs rule).
pub fn clamp_ptrs(p: &mut Pointers, own_len: usize, partner_len: usize) {
    for k in 0..4 {
        if p.pos[k][1] > own_len {
            p.pos[k][1] = own_len;
        }
        if p.pos[k][0] > partner_len {
            p.pos[k][0] = partner_len;
        }
    }
}

fn label_len(s: &Strand, at: usize) -> usize {
    let mut n = 0;
    while is_letter(s.at(at + 1 + n)) {
        n += 1;
    }
    n
}

/// Execute one instruction of the molecule whose pointers are `p` and whose own string is
/// `own`, bound to `partner`. Cleave products are pushed to `products`.
pub fn exec_step(
    p: &mut Pointers,
    own: &mut Strand,
    partner: &mut Strand,
    ctx: &mut Ctx,
    rng: &mut dyn Rng,
    products: &mut Vec<Vec<u8>>,
) -> Outcome {
    macro_rules! s {
        ($side:expr) => {
            if $side == 1 { &*own } else { &*partner }
        };
    }
    let ip = p.get(I);
    let c = s!(p.side[I]).at(ip);
    match c {
        b'$' => {
            let len = label_len(s!(p.side[I]), ip);
            if len == 0 {
                p.side[F] = p.side[I];
                p.set(F, ip);
            } else {
                ctx.label.clear();
                let src = s!(p.side[I]);
                for k in 0..len {
                    ctx.label.push(complement(src.at(ip + 1 + k)));
                }
                let target = s!(p.side[F]).seq();
                let al = ctx.al.align(&ctx.label, target, ctx.tab);
                let prob = search_prob(&al);
                let rno = rng.next_f64() as f32;
                if rno < prob {
                    p.set(F, al.e2);
                } else {
                    p.side[F] = p.side[I];
                    p.set(F, ip + len);
                }
            }
            p.bump(I, 1);
        }
        b'>' => {
            let m = s!(p.side[I]).at(ip + 1);
            let target = match m {
                b'B' => R,
                b'C' => W,
                _ => I,
            };
            let (fs, fp) = (p.side[F], p.get(F));
            p.side[target] = fs;
            p.set(target, fp);
            p.bump(I, 1);
        }
        b'^' => {
            let m = s!(p.side[I]).at(ip + 1);
            let target = match m {
                b'A' => I,
                b'B' => R,
                b'C' => W,
                _ => F,
            };
            p.side[target] ^= 1;
            p.bump(I, 1);
        }
        b'?' => {
            let len = label_len(s!(p.side[I]), ip);
            let rp = p.get(R);
            let at_end = s!(p.side[R]).at(rp) == 0;
            let next = ip + 1;
            let newip = match len {
                0 | 1 => {
                    if at_end { next + len + 1 } else { next + len }
                }
                _ => {
                    ctx.label.clear();
                    ctx.other.clear();
                    let src = s!(p.side[I]);
                    for k in 0..len {
                        ctx.label.push(complement(src.at(next + k)));
                    }
                    let rs = s!(p.side[R]);
                    for k in 0..len {
                        let ch = rs.at(rp + k);
                        if ch == 0 {
                            break;
                        }
                        ctx.other.push(ch);
                    }
                    let al = ctx.al.align(&ctx.label, &ctx.other, ctx.tab);
                    let prob = label_prob(&al, len);
                    let rno = rng.next_f64() as f32;
                    if rno < prob { next + len + 1 } else { next + len }
                }
            };
            p.set(I, newip);
        }
        b'=' => return copy(p, own, partner, ctx, rng),
        b'%' => return cleave(p, own, partner, products, &mut ctx.nascent),
        0 | b'}' => return Outcome::Ended,
        _ => p.bump(I, 1),
    }
    Outcome::Continue
}

fn copy(p: &mut Pointers, own: &mut Strand, partner: &mut Strand, ctx: &mut Ctx, rng: &mut dyn Rng) -> Outcome {
    let maxl = ctx.params.maxl;
    if ctx.nascent.is_some() {
        return copy_nascent(p, own, partner, ctx, rng);
    }
    if p.get(W) >= maxl || p.get(R) >= maxl {
        p.bump(I, 1);
        return Outcome::Ended;
    }
    let ch = if p.side[R] == 1 { own.at(p.get(R)) } else { partner.at(p.get(R)) };
    let write = |p: &Pointers, c: u8, own: &mut Strand, partner: &mut Strand| {
        let at = p.get(W);
        if p.side[W] == 1 { own.put(at, c) } else { partner.put(at, c) }
    };
    if ch != 0 {
        let m = ctx.params.mutate;
        let rno = rng.next_f64() as f32;
        if rno < m {
            let r2 = rng.next_f64() as f32;
            if r2 < 0.5 {
                // insertion: copy, then a random symbol after it
                write(p, ch, own, partner);
                p.bump(W, 1);
                let k = ((rng.next_f64() as f32) * NSYM as f32) as usize;
                write(p, KEY[k.min(NSYM - 1)], own, partner);
                p.bump(W, 1);
            } else if ctx.params.compat_deletion_skip {
                p.bump(I, 1);
            }
            p.bump(R, 1);
        } else {
            let out = if rno < m + m {
                let nb = neighbours(ch, ctx.tab);
                let r3 = rng.next_f64() as f32;
                if (2.0f32 * r3).ceil() as usize <= 1 { nb[0] } else { nb[1] }
            } else {
                ch
            };
            write(p, out, own, partner);
            p.bump(W, 1);
            p.bump(R, 1);
        }
    }
    p.bump(I, 1);
    Outcome::Continue
}

fn cleave(p: &mut Pointers, own: &mut Strand, partner: &mut Strand, products: &mut Vec<Vec<u8>>, nascent: &mut Option<Vec<u8>>) -> Outcome {
    // STRANDCOPY: a cut releases the strand being built; with nothing built it cuts at F as before
    if let Some(n) = nascent.as_mut() {
        if !n.is_empty() {
            products.push(std::mem::take(n));
            p.bump(I, 1);
            return Outcome::Continue;
        }
    }
    let fs = p.side[F];
    let fpos = p.get(F);
    let site = if fs == 1 { &mut *own } else { &mut *partner };
    let len = site.len();
    if fpos < len {
        products.push(site.buf[fpos..len].to_vec());
        for b in &mut site.buf[fpos..len] {
            *b = 0;
        }
        // old rewind_bad_ptrs: pointers past the end go to the end; an empty string's
        // pointers go to its start and everything toggles to the other string
        let plen = partner.len();
        let alen = own.len();
        if plen > 0 {
            for k in 0..4 {
                if p.pos[k][0] > plen {
                    p.pos[k][0] = plen;
                }
            }
        } else {
            for k in 0..4 {
                p.pos[k][0] = 0;
                p.side[k] = 1;
            }
        }
        if alen > 0 {
            for k in 0..4 {
                if p.pos[k][1] > alen {
                    p.pos[k][1] = alen;
                }
            }
        } else {
            for k in 0..4 {
                p.pos[k][1] = 0;
            }
            if plen > 0 {
                p.side = [0; 4];
            }
        }
        if alen == 0 {
            return Outcome::ActiveDestroyed;
        }
        if plen == 0 {
            return Outcome::PassiveDestroyed;
        }
    }
    p.bump(I, 1);
    Outcome::Continue
}

/// STRANDCOPY `=`: append the symbol under R (with the same mutation draws as in-place copying)
/// to the nascent strand. Neither string in the reaction changes.
fn copy_nascent(p: &mut Pointers, own: &mut Strand, partner: &mut Strand, ctx: &mut Ctx, rng: &mut dyn Rng) -> Outcome {
    let maxl = ctx.params.maxl;
    let n = ctx.nascent.as_mut().unwrap();
    if n.len() >= maxl || p.get(R) >= maxl {
        p.bump(I, 1);
        return Outcome::Ended;
    }
    // a template has to be another molecule: with R on the executor's own string, `=` does nothing
    if p.side[R] == 1 {
        p.bump(I, 1);
        return Outcome::Continue;
    }
    let _ = own;
    let mut ch = partner.at(p.get(R));
    let mut step = 1;
    // TRANSLATE: with W toggled off the executor's own string, two letters under R are read as a
    // codon and one symbol is written; an opcode passes through unchanged
    if let (Some(table), 0) = (ctx.params.codon.as_ref(), p.side[W]) {
        let c2 = partner.at(p.get(R) + 1);
        if is_letter(ch) && is_letter(c2) {
            ch = table[(ch - b'A') as usize * 26 + (c2 - b'A') as usize];
            step = 2;
        } else if ctx.params.translate_stop && ch != 0 {
            p.bump(I, 1);
            return Outcome::Ended;
        }
    }
    let n = ctx.nascent.as_mut().unwrap();
    if ch != 0 {
        let m = ctx.params.mutate;
        let rno = rng.next_f64() as f32;
        if rno < m {
            let r2 = rng.next_f64() as f32;
            if r2 < 0.5 {
                n.push(ch);
                let k = ((rng.next_f64() as f32) * NSYM as f32) as usize;
                n.push(KEY[k.min(NSYM - 1)]);
            }
            p.bump(R, step);
        } else {
            let out = if rno < m + m {
                let nb = neighbours(ch, ctx.tab);
                let r3 = rng.next_f64() as f32;
                if (2.0f32 * r3).ceil() as usize <= 1 { nb[0] } else { nb[1] }
            } else {
                ch
            };
            n.push(out);
            p.bump(R, step);
        }
    }
    p.bump(I, 1);
    Outcome::Continue
}
