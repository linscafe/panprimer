//! `pangenome_primer._scan` -- the Phase 2 genome-wide primer matcher.
//!
//! Replaces the `bwa aln` + `bwa samse` pair that `bwa_backend.py` shells out to. That path
//! needs a 5.3 GB on-disk index and 4.45 GB of RAM to answer a query that returns a handful
//! of records; this one streams the 0.9 GB BGZF FASTA that HPRC already distributes and
//! keeps no index at all, which is what makes deleting the `.fa` and the bwa index possible
//! (Phase 3).
//!
//! # The contract
//!
//! Output is defined to be indistinguishable from `binding.find_binding_sites_naive` run
//! over every contig, for the same primers and the same mismatch budget:
//!
//! * `start`/`end` are contig coordinates and `end - start == len(primer)`;
//! * `mismatch_offsets_3p` counts from the primer 3' end (0 = 3' terminal base), which on a
//!   MINUS site is genomic position `start`;
//! * a mismatch is a plain inequality of upper-cased characters -- `N` is never
//!   special-cased, and `N` vs `N` is equal (see `seq.rs`);
//! * no indels;
//! * emission order per primer is contig order, PLUS sites then MINUS sites, ascending by
//!   start within each group.
//!
//! Unlike `bwa aln`, this is exhaustive: no hit cap, no seeding heuristic. See
//! `scanner.rs`.

use std::collections::BTreeSet;

use pyo3::exceptions::{PyIOError, PyValueError};
use pyo3::prelude::*;
use pyo3::types::{PyDict, PyList};

mod bgzf;
mod pool;
mod scanner;
mod seq;

use scanner::{Hit, Index};
use seq::{Contig, FastaParser};

/// (primer index -> its unique sequence), and the encoded forms.
fn prepare(seqs: Vec<String>) -> (Vec<String>, Vec<Vec<u8>>) {
    let uniq: BTreeSet<String> = seqs.into_iter().collect();
    let names: Vec<String> = uniq.into_iter().collect();
    let encoded = names.iter().map(|s| seq::encode(s)).collect();
    (names, encoded)
}

/// Turn per-contig hit lists into `{primer_sequence: [site_tuple, ...]}`.
fn to_python<'py>(
    py: Python<'py>,
    names: &[String],
    index: &Index,
    per_contig: &[(String, Vec<Hit>)],
) -> PyResult<Bound<'py, PyDict>> {
    let out = PyDict::new(py);
    for (pidx, name) in names.iter().enumerate() {
        let sites = PyList::empty(py);
        for (chrom, hits) in per_contig {
            // Two passes so the order is PLUS-then-MINUS within each contig, matching the
            // naive backend's `plus + minus` concatenation exactly.
            for want_minus in [false, true] {
                for hit in hits {
                    let (hp, minus) = index.query_primer(hit.qid);
                    if hp as usize != pidx || minus != want_minus {
                        continue;
                    }
                    let l = index.query_len(hit.qid) as u64;
                    let start = hit.start as u64;
                    sites.append((
                        chrom.as_str(),
                        start,
                        start + l,
                        if minus { "-" } else { "+" },
                        hit.mismatches as u32,
                        index.offsets_3p(hit.qid, hit),
                    ))?;
                }
            }
        }
        out.set_item(name, sites)?;
    }
    Ok(out)
}

/// Genome-wide search over a BGZF (or plain) FASTA.
///
/// `slop` is accepted for signature compatibility with the historical bwa backend, which
/// widened the window fetched around each candidate before re-scoring. This backend scans
/// every position exhaustively, so there is no candidate window to widen and the value is
/// ignored.
///
/// `threads` sizes the crate's rayon pool the first time any scan runs; see `pool.rs`.
#[pyfunction]
#[pyo3(signature = (seqs, path, max_mismatches, slop = 3, threads = None))]
fn scan<'py>(
    py: Python<'py>,
    seqs: Vec<String>,
    path: &str,
    max_mismatches: usize,
    slop: usize,
    threads: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let _ = slop;
    let (names, encoded) = prepare(seqs);
    if names.is_empty() {
        return Ok(PyDict::new(py));
    }
    let index = Index::build(&encoded, max_mismatches).map_err(PyValueError::new_err)?;

    let mut per_contig: Vec<(String, Vec<Hit>)> = Vec::new();
    // `install` makes this the current pool for the whole traversal, so the nested inflate
    // region inside `bgzf::stream` lands on it too rather than on rayon's global pool.
    let res: Result<(), std::io::Error> = py.detach(|| {
        pool::pool(threads).install(|| {
            let mut parser = FastaParser::new(|c: Contig| {
                let hits = index.scan_contig(&c.codes);
                per_contig.push((c.name, hits));
            });
            bgzf::stream(path, |chunk| parser.feed(chunk))?;
            parser.finish();
            Ok(())
        })
    });
    res.map_err(|e| PyIOError::new_err(format!("{path}: {e}")))?;

    to_python(py, &names, &index, &per_contig)
}

/// Same search against an in-memory sequence, for one named contig.
///
/// This is the surface the differential test in
/// `tests/test_rust_backend_differential.py` drives against
/// `binding.find_binding_sites_naive`, and what `binding.find_binding_sites(backend="rust")`
/// dispatches to. Identical code path as `scan` apart from the FASTA/BGZF front end.
#[pyfunction]
#[pyo3(signature = (seqs, ref_seq, max_mismatches, chrom = "chr", threads = None))]
fn scan_seq<'py>(
    py: Python<'py>,
    seqs: Vec<String>,
    ref_seq: &str,
    max_mismatches: usize,
    chrom: &str,
    threads: Option<usize>,
) -> PyResult<Bound<'py, PyDict>> {
    let (names, encoded) = prepare(seqs);
    if names.is_empty() {
        return Ok(PyDict::new(py));
    }
    let index = Index::build(&encoded, max_mismatches).map_err(PyValueError::new_err)?;
    let codes = seq::encode(ref_seq);
    let hits = py.detach(|| pool::pool(threads).install(|| index.scan_contig(&codes)));
    to_python(py, &names, &index, &[(chrom.to_string(), hits)])
}

/// Worker count of the scan pool -- what is live if a scan has run, otherwise what would be
/// chosen. Lets a caller report the real number instead of assuming the one it asked for.
#[pyfunction]
#[pyo3(signature = (threads = None))]
fn pool_threads(threads: Option<usize>) -> usize {
    pool::threads(threads)
}

#[pymodule]
fn _scan(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add("__version__", env!("CARGO_PKG_VERSION"))?;
    m.add_function(wrap_pyfunction!(scan, m)?)?;
    m.add_function(wrap_pyfunction!(scan_seq, m)?)?;
    m.add_function(wrap_pyfunction!(pool_threads, m)?)?;
    Ok(())
}
