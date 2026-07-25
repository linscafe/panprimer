//! Exhaustive <=k-mismatch, gap-free primer search over a contig.
//!
//! # Why this is exact and not a heuristic
//!
//! It replaces `bwa aln`, but it is not an aligner. For every genomic position and every
//! primer (in both orientations) it either proves <=k mismatches or proves >k. There is no
//! seed-extension heuristic, no cap on reported hits, no scoring. The seed index below is
//! purely a *filter*; every candidate it emits is verified base by base with the same rule
//! `binding.find_binding_sites_naive` uses, and the filter is proven never to drop a true
//! site (see `pigeonhole` note).
//!
//! # Pigeonhole
//!
//! Split a primer of length `L` into a prefix half and a suffix half. If the whole primer
//! has <= k mismatches then at least one half has <= k/2 (integer division) -- otherwise
//! both halves would carry >= k/2 + 1 and the total would exceed k. So indexing every
//! variant of the leading `s`-mer and the trailing `s`-mer within `k/2` substitutions, with
//! `s <= L/2` so each seed sits inside its own half, cannot miss a site.
//!
//! For the real case (L=20, k=3): `s = 10`, `k/2 = 1`, so 2 seeds x (1 + 3*10) = 62 codes
//! per primer orientation. A 4^10-entry CSR table over ~2.5k set codes gives a ~0.2%
//! per-position candidate rate -- the 3.1 Gb genome yields a few million verifications,
//! which is nothing next to the ~3.1e9 rolling-hash steps.
//!
//! # Strands
//!
//! A MINUS site is found the same way `binding.find_binding_sites_naive` finds one: by
//! matching the *reverse complement* of the primer against the top strand. Both
//! orientations are separate entries in the same index, so the genome is traversed once.

use rayon::prelude::*;

use crate::seq::{revcomp_codes, two_bit};

/// Largest seed. 4^12 would be a 67 MB CSR offset table; 4^11 is 16.8 MB, which stays
/// friendly to the cache hierarchy while still being far more selective than needed for
/// the 18-27 nt primers this pipeline designs.
const MAX_SEED: usize = 11;
/// Refuse to build an index larger than this many (code, query, offset) triples. Only
/// reachable with an unusually large `max_mismatches`; the fix is a shorter seed.
const MAX_INDEX_ITEMS: usize = 8 << 20;
/// Hard cap on per-site mismatch bookkeeping (`SmallOffsets` is a fixed array).
pub const MAX_MISMATCHES: usize = 16;

/// One primer in one orientation.
struct Query {
    primer_idx: u32,
    minus: bool,
    codes: Vec<u8>,
}

/// A verified binding site, in query-local terms.
#[derive(Clone, Copy)]
pub struct Hit {
    pub qid: u32,
    pub start: u32,
    pub mismatches: u8,
    /// Indices into the top-strand window, ascending. Converted to the 3'-end convention
    /// at emission time (see `Index::offsets_3p`).
    pub window_idx: [u8; MAX_MISMATCHES],
}

pub struct Index {
    queries: Vec<Query>,
    seed: usize,
    /// CSR over 4^seed codes -> slice of `items`.
    csr_off: Vec<u32>,
    items: Vec<(u32, u32)>, // (qid, offset of the seed within the primer)
    max_len: usize,
    max_mismatches: usize,
}

fn enumerate_variants(base: u64, s: usize, budget: usize, out: &mut Vec<u64>) {
    fn rec(cur: u64, pos: usize, s: usize, left: usize, out: &mut Vec<u64>) {
        if pos == s {
            out.push(cur);
            return;
        }
        rec(cur, pos + 1, s, left, out);
        if left > 0 {
            let shift = 2 * (s - 1 - pos);
            let orig = (cur >> shift) & 3;
            for b in 0u64..4 {
                if b != orig {
                    rec((cur & !(3 << shift)) | (b << shift), pos + 1, s, left - 1, out);
                }
            }
        }
    }
    rec(base, 0, s, budget, out);
}

fn seed_code(codes: &[u8]) -> u64 {
    codes.iter().fold(0u64, |acc, &c| (acc << 2) | two_bit(c))
}

impl Index {
    /// `primers` are the unique primer sequences, already encoded (`seq::encode`).
    pub fn build(primers: &[Vec<u8>], max_mismatches: usize) -> Result<Self, String> {
        if primers.is_empty() {
            return Err("no primer sequences given".into());
        }
        if max_mismatches > MAX_MISMATCHES {
            return Err(format!(
                "max_mismatches={max_mismatches} exceeds the supported limit of {MAX_MISMATCHES}"
            ));
        }
        let mut queries: Vec<Query> = Vec::with_capacity(primers.len() * 2);
        for (i, p) in primers.iter().enumerate() {
            if p.is_empty() {
                return Err("empty primer sequence".into());
            }
            if p.len() > 255 {
                return Err(format!("primer of length {} exceeds 255 nt", p.len()));
            }
            queries.push(Query { primer_idx: i as u32, minus: false, codes: p.clone() });
            queries.push(Query { primer_idx: i as u32, minus: true, codes: revcomp_codes(p) });
        }
        let min_len = queries.iter().map(|q| q.codes.len()).min().unwrap();
        let max_len = queries.iter().map(|q| q.codes.len()).max().unwrap();
        let budget = max_mismatches / 2;

        let mut seed = (min_len / 2).min(MAX_SEED).max(1);
        let (csr_off, items) = loop {
            match Self::try_build_index(&queries, seed, budget) {
                Some(v) => break v,
                None => {
                    if seed == 1 {
                        return Err(format!(
                            "cannot build a seed index for max_mismatches={max_mismatches} \
                             (variant explosion); lower max_mismatches"
                        ));
                    }
                    seed -= 1;
                }
            }
        };
        Ok(Self {
            queries,
            seed,
            csr_off,
            items,
            max_len,
            max_mismatches,
        })
    }

    fn try_build_index(
        queries: &[Query],
        seed: usize,
        budget: usize,
    ) -> Option<(Vec<u32>, Vec<(u32, u32)>)> {
        let table = 1usize << (2 * seed);
        let mut triples: Vec<(u64, u32, u32)> = Vec::new();
        let mut vbuf: Vec<u64> = Vec::new();
        for (qid, q) in queries.iter().enumerate() {
            let l = q.codes.len();
            // Both seeds live inside their own half because seed <= min_len/2 <= l/2.
            let mut offsets = vec![0usize, l - seed];
            offsets.dedup();
            for &off in &offsets {
                vbuf.clear();
                enumerate_variants(seed_code(&q.codes[off..off + seed]), seed, budget, &mut vbuf);
                if triples.len() + vbuf.len() > MAX_INDEX_ITEMS {
                    return None;
                }
                for &code in &vbuf {
                    triples.push((code, qid as u32, off as u32));
                }
            }
        }
        triples.sort_unstable();
        triples.dedup();
        let mut csr_off = vec![0u32; table + 1];
        for &(code, _, _) in &triples {
            csr_off[code as usize + 1] += 1;
        }
        for i in 0..table {
            csr_off[i + 1] += csr_off[i];
        }
        let items: Vec<(u32, u32)> = triples.iter().map(|&(_, q, o)| (q, o)).collect();
        Some((csr_off, items))
    }

    #[inline]
    pub fn query_primer(&self, qid: u32) -> (u32, bool) {
        let q = &self.queries[qid as usize];
        (q.primer_idx, q.minus)
    }

    #[inline]
    pub fn query_len(&self, qid: u32) -> usize {
        self.queries[qid as usize].codes.len()
    }

    /// Convert a hit's top-strand window indices to `mismatch_offsets_3p`, ascending.
    ///
    /// PLUS: window index `k` is primer index `k`, whose distance from the 3' end is
    /// `L-1-k`. Ascending `k` therefore yields descending offsets -- hence the reverse.
    ///
    /// MINUS: `binding.find_binding_sites_naive` scans with the reverse-complemented
    /// primer (offset `L-1-k`) and then relabels every offset to `(L-1) - (L-1-k) = k`.
    /// So for a MINUS site the 3'-offset IS the top-strand window index, ascending
    /// already: offset 0 sits at genomic position `start`.
    pub fn offsets_3p(&self, qid: u32, hit: &Hit) -> Vec<i64> {
        let l = self.query_len(qid);
        let n = hit.mismatches as usize;
        let idx = &hit.window_idx[..n];
        if self.queries[qid as usize].minus {
            idx.iter().map(|&k| k as i64).collect()
        } else {
            idx.iter().rev().map(|&k| (l - 1 - k as usize) as i64).collect()
        }
    }

    /// Scan `[lo, hi)` **of primer start positions** in `codes`, appending verified hits.
    fn scan_range(&self, codes: &[u8], lo: usize, hi: usize, out: &mut Vec<Hit>) {
        let n = codes.len();
        let s = self.seed;
        if n < s || lo >= hi {
            return;
        }
        let mask = if 2 * s >= 64 { u64::MAX } else { (1u64 << (2 * s)) - 1 };
        // A hit at `start` is discovered from the seed window at `start + off`, with
        // `off` in [0, L-s]. To cover every start in [lo, hi) we must therefore visit seed
        // windows in [lo, hi + max_len - s].
        let p_end = (hi + self.max_len.saturating_sub(s)).min(n - s + 1);
        if lo >= p_end {
            return;
        }
        let mut code = 0u64;
        for i in lo..(lo + s - 1).min(n) {
            code = ((code << 2) | two_bit(codes[i])) & mask;
        }
        for p in lo..p_end {
            code = ((code << 2) | two_bit(codes[p + s - 1])) & mask;
            let a = self.csr_off[code as usize] as usize;
            let b = self.csr_off[code as usize + 1] as usize;
            if a == b {
                continue;
            }
            for &(qid, off) in &self.items[a..b] {
                let off = off as usize;
                if off > p {
                    continue;
                }
                let start = p - off;
                if start < lo || start >= hi {
                    continue;
                }
                let q = &self.queries[qid as usize];
                let l = q.codes.len();
                if start + l > n {
                    continue;
                }
                let window = &codes[start..start + l];
                let mut mm = 0usize;
                let mut widx = [0u8; MAX_MISMATCHES];
                let mut over = false;
                for k in 0..l {
                    if window[k] != q.codes[k] {
                        if mm == self.max_mismatches {
                            over = true;
                            break;
                        }
                        widx[mm] = k as u8;
                        mm += 1;
                    }
                }
                if !over {
                    out.push(Hit {
                        qid,
                        start: start as u32,
                        mismatches: mm as u8,
                        window_idx: widx,
                    });
                }
            }
        }
    }

    /// Scan a whole contig. Returned hits are deduplicated on `(qid, start)` -- the two
    /// seeds of one primer can both fire on the same site -- and sorted by `(qid, start)`.
    pub fn scan_contig(&self, codes: &[u8]) -> Vec<Hit> {
        let n = codes.len();
        if n < self.seed {
            return Vec::new();
        }
        const CHUNK: usize = 4 << 20;
        let n_chunks = n.div_ceil(CHUNK).max(1);
        let mut hits: Vec<Hit> = (0..n_chunks)
            .into_par_iter()
            .map(|ci| {
                let lo = ci * CHUNK;
                let hi = ((ci + 1) * CHUNK).min(n);
                let mut v = Vec::new();
                self.scan_range(codes, lo, hi, &mut v);
                v
            })
            .reduce(Vec::new, |mut a, mut b| {
                a.append(&mut b);
                a
            });
        hits.sort_unstable_by_key(|h| (h.qid, h.start));
        hits.dedup_by_key(|h| (h.qid, h.start));
        hits
    }
}
