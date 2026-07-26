//! Minimal streaming BGZF (and plain-FASTA) reader.
//!
//! BGZF is gzip with one extra rule: every member carries an `BC` FEXTRA subfield holding
//! `BSIZE-1`, the total on-disk size of that member. That makes the file splittable without
//! inflating it, which is the whole point here -- we cut a read buffer into whole members and
//! inflate them across a rayon pool, so the ~0.9 GB compressed haplotype is decompressed at
//! aggregate core speed rather than one-core zlib speed.
//!
//! Two deliberate choices:
//!
//! * pure-Rust inflate (`flate2` default features => `miniz_oxide`). Linking htslib or
//!   libdeflate would drag a C toolchain into every source build, and the whole point of the
//!   `maturin` packaging is that a user with a Rust toolchain can build this from source.
//! * bounded streaming rather than mmap. Peak RSS is the metric this phase is judged on
//!   (4.45 GB bwa baseline); an mmap of the whole file would charge every touched page to
//!   RSS and make the measurement meaningless.
//!
//! A file that does not start with the gzip magic is streamed verbatim, so a plain
//! uncompressed `.fa` also works (slower -- 3.1 GB of I/O instead of 0.9 GB).

use std::fs::File;
use std::io::{self, BufReader, Read};

use flate2::{Decompress, FlushDecompress};
use rayon::prelude::*;

/// Bytes of compressed input pulled per batch. Each batch is split into whole BGZF members
/// and inflated in parallel, so this also sets the parallel granularity: 32 MB is ~500
/// members of 64 KB, plenty to keep 16 cores busy, and bounds the transient buffers.
const READ_CHUNK: usize = 32 << 20;

/// Parsed BGZF member header: (total member size on disk, XLEN).
fn parse_member_header(b: &[u8]) -> Option<(usize, usize)> {
    if b.len() < 12 {
        return None;
    }
    if b[0] != 0x1f || b[1] != 0x8b || b[2] != 8 {
        return None;
    }
    if b[3] & 0x04 == 0 {
        return None; // no FEXTRA => not BGZF
    }
    let xlen = u16::from_le_bytes([b[10], b[11]]) as usize;
    if b.len() < 12 + xlen {
        return None;
    }
    let extra = &b[12..12 + xlen];
    let mut i = 0usize;
    while i + 4 <= extra.len() {
        let slen = u16::from_le_bytes([extra[i + 2], extra[i + 3]]) as usize;
        if extra[i] == b'B' && extra[i + 1] == b'C' && slen == 2 && i + 6 <= extra.len() {
            let bsize = u16::from_le_bytes([extra[i + 4], extra[i + 5]]) as usize + 1;
            return Some((bsize, xlen));
        }
        i += 4 + slen;
    }
    None
}

/// Allocate one decompressor. **`#[inline(never)] is load-bearing, not a hint.**
///
/// `Decompress::new` bottoms out in `InflateState::new_boxed`, which is `Box::default()` --
/// it materialises the whole ~44 KB `InflateState` (a 32 KB LZ dictionary plus three Huffman
/// tables) as a stack temporary and then memcpy's it into the box. LLVM does not elide that
/// temporary.
///
/// Inlined into the closure below, that temporary lands in the frame of rayon's
/// `bridge_producer_consumer::helper` -- which is *recursive*, splitting the member list and
/// re-entering itself through `join`. Measured frame was 0xd298 = 53,912 bytes, against a
/// 2 MiB default worker stack: roughly 38 frames to overflow, and a work-stealing worker
/// nests deeper than the split depth alone suggests. That was ISSUE-001, an intermittent
/// SIGSEGV whose frequency depended on the steal schedule rather than on the input.
///
/// Keeping this out of line confines the temporary to a leaf frame that pops immediately;
/// only the boxed state (on the heap) crosses back. See `docs/scanner_notes.md`.
#[inline(never)]
fn new_decompressor() -> Decompress {
    Decompress::new(false) // false => raw deflate, no zlib wrapper
}

/// Inflate one BGZF member with a **borrowed, reused** decompressor.
///
/// `d` is supplied per rayon job by `map_init` and reset per member, so a batch of ~500
/// members costs one 44 KB allocation per worker instead of 500. `reset(false)` restores the
/// raw-deflate state exactly as `new` left it, and also zeroes `total_out`, which the ISIZE
/// check below reads.
///
/// The reset is load-bearing, and was checked by deleting it: the very next member in a job
/// then inflates wrong (measured: 65,276 bytes against an ISIZE of 65,270). The length check
/// below catches that and fails the scan loudly, which is the point of validating ISIZE per
/// member rather than trusting the decompressor.
fn inflate_raw(d: &mut Decompress, src: &[u8], isize_hint: usize) -> io::Result<Vec<u8>> {
    let mut out = vec![0u8; isize_hint];
    if isize_hint == 0 {
        return Ok(out);
    }
    d.reset(false);
    d.decompress(src, &mut out, FlushDecompress::Finish)
        .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, format!("bgzf inflate: {e}")))?;
    let n = d.total_out() as usize;
    if n != isize_hint {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("bgzf member inflated to {n} bytes, ISIZE said {isize_hint}"),
        ));
    }
    Ok(out)
}

/// Stream `path`, handing decompressed byte runs to `sink` **in file order**.
///
/// `sink` sees the exact concatenation of the file's uncompressed content, split at
/// arbitrary boundaries; it must be able to resume mid-line (`FastaParser` is).
pub fn stream<F>(path: &str, mut sink: F) -> io::Result<()>
where
    F: FnMut(&[u8]),
{
    let mut file = BufReader::with_capacity(1 << 20, File::open(path)?);

    // Sniff the magic: a plain (non-gzip) FASTA is passed straight through.
    let mut magic = [0u8; 2];
    let mut got = 0usize;
    while got < 2 {
        let n = file.read(&mut magic[got..])?;
        if n == 0 {
            break;
        }
        got += n;
    }
    if got < 2 || magic[0] != 0x1f || magic[1] != 0x8b {
        sink(&magic[..got]);
        let mut buf = vec![0u8; READ_CHUNK];
        loop {
            let n = file.read(&mut buf)?;
            if n == 0 {
                return Ok(());
            }
            sink(&buf[..n]);
        }
    }

    let mut buf: Vec<u8> = Vec::with_capacity(READ_CHUNK + (1 << 17));
    buf.extend_from_slice(&magic);
    let mut eof = false;
    loop {
        if !eof {
            let start = buf.len();
            buf.resize(start + READ_CHUNK, 0);
            let mut n = 0usize;
            while n < READ_CHUNK {
                let r = file.read(&mut buf[start + n..])?;
                if r == 0 {
                    eof = true;
                    break;
                }
                n += r;
            }
            buf.truncate(start + n);
        }
        if buf.is_empty() {
            return Ok(());
        }

        // Cut the buffer into whole members: (cdata_start, cdata_end, ISIZE).
        let mut members: Vec<(usize, usize, usize)> = Vec::new();
        let mut pos = 0usize;
        while pos < buf.len() {
            let (total, xlen) = match parse_member_header(&buf[pos..]) {
                Some(v) => v,
                None => {
                    if buf.len() - pos >= 64 {
                        return Err(io::Error::new(
                            io::ErrorKind::InvalidData,
                            format!(
                                "{path}: not a BGZF file (no BC subfield at byte {pos}). \
                                 HPRC .fa.gz files are BGZF; a plain `gzip` file is not \
                                 seekable/splittable and is not supported."
                            ),
                        ));
                    }
                    break; // truncated header: carry it to the next read
                }
            };
            if pos + total > buf.len() {
                break; // truncated member: carry
            }
            let isize_ = u32::from_le_bytes(
                buf[pos + total - 4..pos + total].try_into().unwrap(),
            ) as usize;
            members.push((pos + 12 + xlen, pos + total - 8, isize_));
            pos += total;
        }

        if members.is_empty() {
            if eof {
                if pos < buf.len() {
                    return Err(io::Error::new(
                        io::ErrorKind::UnexpectedEof,
                        format!("{path}: truncated BGZF member at end of file"),
                    ));
                }
                return Ok(());
            }
            continue; // need more bytes before a whole member is available
        }

        let inflated: Vec<Vec<u8>> = members
            .par_iter()
            .map_init(new_decompressor, |d, &(cs, ce, isz)| {
                inflate_raw(d, &buf[cs..ce], isz)
            })
            .collect::<io::Result<Vec<_>>>()?;
        for block in &inflated {
            if !block.is_empty() {
                sink(block);
            }
        }
        drop(inflated);

        buf.drain(..pos);
        if eof && buf.is_empty() {
            return Ok(());
        }
        if eof && pos == 0 {
            return Err(io::Error::new(
                io::ErrorKind::UnexpectedEof,
                format!("{path}: truncated BGZF member at end of file"),
            ));
        }
    }
}
