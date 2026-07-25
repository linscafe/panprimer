//! Base encoding and streaming FASTA parsing.
//!
//! # The encoding, and why it is what it is
//!
//! `binding.find_binding_sites_naive` decides "is this a mismatch?" with a plain Python
//! `!=` between two **upper-cased characters**. That has three consequences this backend
//! must reproduce exactly, because `tests/test_binding_naive.py` pins all three:
//!
//! * soft-masked lowercase reference bases are equal to their uppercase primer base;
//! * `N` in the reference is a mismatch against any real base (never special-cased);
//! * `N` in the primer against `N` in the reference is **equal**, i.e. NOT a mismatch
//!   (`test_n_vs_n_is_literal_equality_not_a_mismatch`). The same goes for any other
//!   IUPAC ambiguity code compared against itself -- and `R` vs `N` IS a mismatch.
//!
//! So `code()` maps A/C/G/T (either case) to 0..=3 and leaves every other byte as its
//! upper-cased self. Since no non-ACGT byte in a FASTA can be < 4, equality of codes is
//! exactly equality of upper-cased characters. That single rule reproduces all three
//! behaviours without any of them being a special case here either.

/// A/C/G/T (either case) -> 0..=3; anything else -> its upper-cased byte (always >= 32).
pub const CODE: [u8; 256] = {
    let mut t = [0u8; 256];
    let mut i = 0usize;
    while i < 256 {
        let b = i as u8;
        t[i] = match b {
            b'A' | b'a' => 0,
            b'C' | b'c' => 1,
            b'G' | b'g' => 2,
            b'T' | b't' => 3,
            _ => b.to_ascii_uppercase(),
        };
        i += 1;
    }
    t
};

/// 2-bit projection used for *seeding only*. Ambiguous bases collapse to 0 (`A`).
///
/// This never causes a missed site. Equal characters always have equal codes and therefore
/// equal 2-bit projections, so the number of differing 2-bit symbols in a window is <= the
/// number of true character mismatches. The seed index is built from every variant within
/// `max_mismatches/2` differences, so a window that passes the pigeonhole bound on true
/// mismatches also passes it on 2-bit symbols. It only ever over-generates candidates,
/// which the exact verifier then rejects.
#[inline(always)]
pub fn two_bit(code: u8) -> u64 {
    if code < 4 {
        code as u64
    } else {
        0
    }
}

pub fn encode(seq: &str) -> Vec<u8> {
    seq.as_bytes().iter().map(|&b| CODE[b as usize]).collect()
}

/// Reverse complement of an already-encoded sequence, in code space.
/// Complementing is `3 - c` for A/C/G/T; anything ambiguous is complemented as
/// `str.maketrans("ACGTNacgtn", "TGCANtgcan")` does -- i.e. `N`->`N`, everything else
/// left alone -- matching `binding.revcomp`.
pub fn revcomp_codes(codes: &[u8]) -> Vec<u8> {
    codes
        .iter()
        .rev()
        .map(|&c| if c < 4 { 3 - c } else { c })
        .collect()
}

/// One fully-decoded contig: name as pysam/`.fai` report it (header up to first
/// whitespace) and the encoded top-strand sequence.
pub struct Contig {
    pub name: String,
    pub codes: Vec<u8>,
}

/// Incremental FASTA parser. `feed` may be called with arbitrary byte splits (a BGZF member
/// boundary lands mid-line constantly); `finish` closes the last contig.
///
/// Each completed contig is handed to the callback the moment its `>` terminator is seen,
/// so the scanner can consume it and free it while decompression continues -- peak RSS is
/// one contig (251 MB for the largest HPRC chromosome) rather than a whole genome.
pub struct FastaParser<F: FnMut(Contig)> {
    on_contig: F,
    in_header: bool,
    at_line_start: bool,
    header: Vec<u8>,
    name: Option<String>,
    codes: Vec<u8>,
}

impl<F: FnMut(Contig)> FastaParser<F> {
    pub fn new(on_contig: F) -> Self {
        Self {
            on_contig,
            in_header: false,
            at_line_start: true,
            header: Vec::new(),
            name: None,
            codes: Vec::new(),
        }
    }

    fn flush(&mut self) {
        if let Some(name) = self.name.take() {
            let codes = std::mem::take(&mut self.codes);
            (self.on_contig)(Contig { name, codes });
        }
    }

    fn open_contig(&mut self) {
        // pysam / samtools faidx name a reference by the header up to the first whitespace.
        let h = String::from_utf8_lossy(&self.header);
        let name = h.split_whitespace().next().unwrap_or("").to_string();
        self.name = Some(name);
        self.header.clear();
    }

    pub fn feed(&mut self, mut chunk: &[u8]) {
        while !chunk.is_empty() {
            if self.in_header {
                match memchr::memchr(b'\n', chunk) {
                    Some(i) => {
                        self.header.extend_from_slice(&chunk[..i]);
                        if self.header.last() == Some(&b'\r') {
                            self.header.pop();
                        }
                        self.in_header = false;
                        self.at_line_start = true;
                        self.open_contig();
                        chunk = &chunk[i + 1..];
                    }
                    None => {
                        self.header.extend_from_slice(chunk);
                        return;
                    }
                }
                continue;
            }

            if self.at_line_start && chunk[0] == b'>' {
                self.flush();
                self.in_header = true;
                self.header.clear();
                chunk = &chunk[1..];
                continue;
            }

            // Bulk path: everything up to the next newline is sequence.
            let end = memchr::memchr(b'\n', chunk).unwrap_or(chunk.len());
            let mut run = &chunk[..end];
            if run.last() == Some(&b'\r') {
                run = &run[..run.len() - 1];
            }
            if self.name.is_some() && !run.is_empty() {
                self.codes.reserve(run.len());
                self.codes.extend(run.iter().map(|&b| CODE[b as usize]));
            }
            if end < chunk.len() {
                self.at_line_start = true;
                chunk = &chunk[end + 1..];
            } else {
                self.at_line_start = false;
                return;
            }
        }
    }

    pub fn finish(mut self) {
        if self.in_header {
            self.open_contig();
        }
        self.flush();
    }
}
