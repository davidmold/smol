//! Config files, in the old Stringmol format: `KEY value` lines, `%` or `#` comments,
//! `AGENT seq [count] ...` followed by `GRIDPOS x y` lines.

use std::collections::HashMap;

#[derive(Clone, Debug)]
pub struct Config {
    pub gridx: usize,
    pub gridy: usize,
    pub decay: f64,
    pub mutate: f32,
    pub maxl: usize,
    pub seed: u64,
    pub nsteps: u64,
    pub popdy_every: u64,
    pub snap_every: u64,
    pub copy_every: u64,
    pub ckpt_every: u64,
    pub keep_ckpt: usize,
    pub bindlen: usize,
    /// > 0: absolute-score binding rule with this scale (see align::bind_prob_rule)
    pub bind_t: f64,
    /// Old partner choice: always the first free neighbour in scan order (quirk 1).
    pub compat_partner: bool,
    /// Decay in pairs. Default (as the old code): one roll per pair, both die.
    /// DECAY_INDEPENDENT 1: each member rolls; a decaying molecule dies, its partner is freed.
    pub decay_independent: bool,
    /// A failed encounter also uses up the partner's turn (old rule; default 0 = it does not).
    pub encounter_turn: bool,
    /// Strand displacement (design note): a free molecule may break into a reaction.
    pub displace: bool,
    /// Copy instructions after which a reaction can no longer be broken into.
    pub displace_k: f64,
    /// Overall scale on the displacement probability.
    pub displace_rate: f64,
    /// On displacement, undo the abandoned reaction's writes to both strands.
    pub displace_restore: bool,
    /// Two-sided reactions: the passive also runs its own code, one instruction after the
    /// active's each step.
    pub twosided: bool,
    /// Copies are built on a separate nascent strand (see reaction::copy_nascent).
    pub strandcopy: bool,
    /// (sequence, x, y)
    pub agents: Vec<(Vec<u8>, usize, usize)>,
    /// Every key/value as read, for run.json.
    pub raw: Vec<(String, String)>,
}

fn strip(line: &str) -> &str {
    let l = line.trim();
    if l.starts_with('%') || l.starts_with('#') {
        ""
    } else {
        l
    }
}

impl Config {
    pub fn parse(text: &str) -> Result<Config, String> {
        let mut kv: HashMap<String, String> = HashMap::new();
        let mut raw = Vec::new();
        let mut agents = Vec::new();
        let mut cur: Option<Vec<u8>> = None;
        for line in text.lines() {
            let l = strip(line);
            if l.is_empty() {
                continue;
            }
            let mut it = l.split_whitespace();
            let key = it.next().unwrap().to_string();
            let rest: Vec<&str> = it.collect();
            match key.as_str() {
                "AGENT" => {
                    let seq = rest.first().ok_or("AGENT without a sequence")?;
                    cur = Some(seq.as_bytes().to_vec());
                }
                "GRIDPOS" => {
                    let seq = cur.clone().ok_or("GRIDPOS before any AGENT")?;
                    let x: usize = rest.first().and_then(|v| v.parse().ok()).ok_or("bad GRIDPOS")?;
                    let y: usize = rest.get(1).and_then(|v| v.parse().ok()).ok_or("bad GRIDPOS")?;
                    agents.push((seq, x, y));
                }
                _ => {
                    if let Some(v) = rest.first() {
                        kv.insert(key.clone(), v.to_string());
                        raw.push((key, v.to_string()));
                    }
                }
            }
        }
        let get = |k: &str| kv.get(k).cloned();
        let num = |k: &str, d: f64| -> Result<f64, String> {
            match get(k) {
                Some(v) => v.parse::<f64>().map_err(|_| format!("bad value for {k}: {v}")),
                None => Ok(d),
            }
        };
        let mutate = match get("MUTATE") {
            Some(v) => v.parse::<f32>().map_err(|_| format!("bad MUTATE {v}"))?,
            None => return Err("MUTATE is required (spec: Copying and mutation)".into()),
        };
        let gridx = num("GRIDX", 125.0)? as usize;
        let gridy = num("GRIDY", 100.0)? as usize;
        for (_, x, y) in &agents {
            if *x >= gridx || *y >= gridy {
                return Err(format!("GRIDPOS {x} {y} is outside the {gridx}x{gridy} grid"));
            }
        }
        Ok(Config {
            gridx,
            gridy,
            decay: num("DECAY", 0.0005)?,
            mutate,
            maxl: num("MAXLEN", 2000.0)? as usize,
            seed: num("RANDSEED", 1.0)? as u64,
            nsteps: num("NSTEPS", 1_000_000.0)? as u64,
            popdy_every: num("POPDYEVERY", 100.0)? as u64,
            snap_every: num("SNAPEVERY", 500.0)? as u64,
            copy_every: num("COPYEVERY", 1000.0)? as u64,
            ckpt_every: num("CKPTEVERY", 50_000.0)? as u64,
            keep_ckpt: num("KEEPCKPT", 3.0)? as usize,
            bindlen: num("BINDLEN", 0.0)? as usize,
            bind_t: num("BIND_T", 0.0)?,
            compat_partner: num("COMPAT_PARTNER", 0.0)? != 0.0,
            decay_independent: num("DECAY_INDEPENDENT", 0.0)? != 0.0,
            encounter_turn: num("ENCOUNTER_TURN", 0.0)? != 0.0,
            displace: num("DISPLACE", 0.0)? != 0.0,
            displace_k: num("DISPLACE_K", 10.0)?,
            displace_rate: num("DISPLACE_RATE", 1.0)?,
            displace_restore: num("DISPLACE_RESTORE", 0.0)? != 0.0,
            twosided: num("TWOSIDED", 0.0)? != 0.0,
            strandcopy: num("STRANDCOPY", 0.0)? != 0.0,
            agents,
            raw,
        })
    }
}
