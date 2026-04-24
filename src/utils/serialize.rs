/// Binary serialization for commitment and proof files.
///
/// Both formats are little-endian (native x86-64).  All `[u32]` payloads are
/// written/read via [`bytemuck::cast_slice`] for zero-copy conversion.
use std::fs::File;
use std::io::{self, BufReader, BufWriter, Read, Write};
use std::path::Path;

use bytemuck::cast_slice;

use crate::hachi::commit::Commitment;
use crate::hachi::prove::SumcheckProof;
use crate::utils::ds::{Align64, AlignedU32Vec, AlignedU8Vec};

// ── helpers ──────────────────────────────────────────────────────────────────

fn write_u32_slice(w: &mut impl Write, data: &[u32]) -> io::Result<()> {
    w.write_all(cast_slice::<u32, u8>(data))
}

fn read_u32_vec(r: &mut impl Read, count: usize) -> io::Result<Vec<u32>> {
    let mut buf = vec![0u32; count];
    r.read_exact(bytemuck::cast_slice_mut::<u32, u8>(&mut buf))?;
    Ok(buf)
}

fn write_u32(w: &mut impl Write, v: u32) -> io::Result<()> {
    w.write_all(&v.to_le_bytes())
}

fn read_u32(r: &mut impl Read) -> io::Result<u32> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    Ok(u32::from_le_bytes(buf))
}

fn write_magic(w: &mut impl Write, magic: &[u8; 4]) -> io::Result<()> {
    w.write_all(magic)
}

fn read_magic(r: &mut impl Read, expected: &[u8; 4]) -> io::Result<()> {
    let mut buf = [0u8; 4];
    r.read_exact(&mut buf)?;
    if &buf != expected {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("bad magic: expected {:?}, got {:?}", expected, buf),
        ));
    }
    Ok(())
}

// ── AlignedU32Vec helpers ────────────────────────────────────────────────────

/// Build an [`AlignedU32Vec`] from a flat `[u32]` slice.
fn aligned_u32_vec_from_slice(data: &[u32]) -> AlignedU32Vec {
    let len = data.len();
    let chunks: Vec<Align64<[u32; 16]>> = data
        .chunks(16)
        .map(|c| {
            let mut arr = [0u32; 16];
            arr[..c.len()].copy_from_slice(c);
            Align64(arr)
        })
        .collect();
    AlignedU32Vec { inner: chunks, len }
}

// ── Commitment I/O ──────────────────────────────────────────────────────────

const COMMITMENT_MAGIC: &[u8; 4] = b"HCMT";
const COMMITMENT_VERSION: u32 = 1;

pub fn write_commitment(path: &Path, commitment: &Commitment) -> io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);

    write_magic(&mut w, COMMITMENT_MAGIC)?;
    write_u32(&mut w, COMMITMENT_VERSION)?;
    write_u32(&mut w, commitment.u.len as u32)?;
    write_u32(&mut w, commitment.r.len as u32)?;
    write_u32(&mut w, commitment.t.len as u32)?;

    write_u32_slice(&mut w, &commitment.u)?;
    write_u32_slice(&mut w, &commitment.r)?;
    write_u32_slice(&mut w, &commitment.t)?;

    w.flush()
}

pub fn read_commitment(path: &Path) -> io::Result<Commitment> {
    let f = File::open(path)?;
    let mut r = BufReader::new(f);

    read_magic(&mut r, COMMITMENT_MAGIC)?;
    let version = read_u32(&mut r)?;
    if version != COMMITMENT_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported commitment version: {version}"),
        ));
    }

    let u_len = read_u32(&mut r)? as usize;
    let r_len = read_u32(&mut r)? as usize;
    let t_len = read_u32(&mut r)? as usize;

    let u = aligned_u32_vec_from_slice(&read_u32_vec(&mut r, u_len)?);
    let rv = aligned_u32_vec_from_slice(&read_u32_vec(&mut r, r_len)?);
    let t = aligned_u32_vec_from_slice(&read_u32_vec(&mut r, t_len)?);

    Ok(Commitment { u, r: rv, t })
}

// ── Proof I/O ───────────────────────────────────────────────────────────────

const PROOF_MAGIC: &[u8; 4] = b"HPRF";
const PROOF_VERSION: u32 = 1;

/// Bundle returned by [`read_proof`].
pub struct ProofBundle {
    pub proof: SumcheckProof,
    pub target_sum: [u32; 4],
}

/// Write a [`SumcheckProof`] together with the target sum.
pub fn write_proof(
    path: &Path,
    proof: &SumcheckProof,
    target_sum: [u32; 4],
) -> io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);

    write_magic(&mut w, PROOF_MAGIC)?;
    write_u32(&mut w, PROOF_VERSION)?;
    write_u32(&mut w, proof.p_alg_x.len() as u32)?;
    write_u32(&mut w, proof.p_alg_y.len() as u32)?;

    // X rounds
    write_vec_arr::<3>(&mut w, &proof.p_alg_x)?;
    write_vec_arr::<7>(&mut w, &proof.p_norm_x)?;
    write_vec_field(&mut w, &proof.r_x)?;

    // Y rounds
    write_vec_arr::<3>(&mut w, &proof.p_alg_y)?;
    write_vec_arr::<7>(&mut w, &proof.p_norm_y)?;
    write_vec_field(&mut w, &proof.r_y)?;

    // Final evaluation
    write_u32_slice(&mut w, &proof.final_w)?;

    // Target sum
    write_u32_slice(&mut w, &target_sum)?;

    w.flush()
}

pub fn read_proof(path: &Path) -> io::Result<ProofBundle> {
    let f = File::open(path)?;
    let mut r = BufReader::new(f);

    read_magic(&mut r, PROOF_MAGIC)?;
    let version = read_u32(&mut r)?;
    if version != PROOF_VERSION {
        return Err(io::Error::new(
            io::ErrorKind::InvalidData,
            format!("unsupported proof version: {version}"),
        ));
    }

    let rounds_x = read_u32(&mut r)? as usize;
    let rounds_y = read_u32(&mut r)? as usize;

    // X rounds
    let p_alg_x = read_vec_arr::<3>(&mut r, rounds_x)?;
    let p_norm_x = read_vec_arr::<7>(&mut r, rounds_x)?;
    let r_x = read_vec_field(&mut r, rounds_x)?;

    // Y rounds
    let p_alg_y = read_vec_arr::<3>(&mut r, rounds_y)?;
    let p_norm_y = read_vec_arr::<7>(&mut r, rounds_y)?;
    let r_y = read_vec_field(&mut r, rounds_y)?;

    // Final evaluation
    let final_w = read_field4(&mut r)?;

    let proof = SumcheckProof {
        p_alg_x,
        p_norm_x,
        r_x,
        p_alg_y,
        p_norm_y,
        r_y,
        final_w,
    };

    // Target sum
    let target_sum = read_field4(&mut r)?;

    Ok(ProofBundle {
        proof,
        target_sum,
    })
}

// ── Witness I/O ─────────────────────────────────────────────────────────────

/// Write an [`AlignedU8Vec`] as a raw binary witness file.
///
/// Only the logical `len` prefix is written — any trailing padding in the
/// final 64-byte chunk is discarded.
pub fn write_witness(path: &Path, witness: &AlignedU8Vec) -> io::Result<()> {
    let f = File::create(path)?;
    let mut w = BufWriter::new(f);
    w.write_all(&witness[..witness.len])?;
    w.flush()
}

/// Read a raw binary witness file into an [`AlignedU8Vec`].
///
/// Allocates the 64-byte-aligned backing store once (zero-initialised, so any
/// tail past `len` is padding) and reads the file straight into it — no
/// intermediate `Vec<u8>` and no per-chunk copy.
pub fn read_witness(path: &Path) -> io::Result<AlignedU8Vec> {
    let mut f = File::open(path)?;
    let len = f.metadata()?.len() as usize;
    let num_chunks = len.div_ceil(64);
    let mut inner: Vec<Align64<[u8; 64]>> = vec![Align64([0u8; 64]); num_chunks];

    // SAFETY: `Align64<[u8; 64]>` is `#[repr(C, align(64))]` over `[u8; 64]`,
    // so `num_chunks * 64` bytes of `inner` are valid, initialised storage.
    // We read `len <= num_chunks * 64` bytes; the remainder stays zero.
    let buf: &mut [u8] = unsafe {
        std::slice::from_raw_parts_mut(inner.as_mut_ptr() as *mut u8, len)
    };
    f.read_exact(buf)?;

    Ok(AlignedU8Vec { inner, len })
}

// ── Vec<[[u32; 4]; N]> helpers ──────────────────────────────────────────────

fn write_vec_arr<const N: usize>(
    w: &mut impl Write,
    data: &[[[u32; 4]; N]],
) -> io::Result<()> {
    for round in data {
        for elem in round {
            write_u32_slice(w, elem)?;
        }
    }
    Ok(())
}

fn read_vec_arr<const N: usize>(
    r: &mut impl Read,
    count: usize,
) -> io::Result<Vec<[[u32; 4]; N]>> {
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        let mut round = [[0u32; 4]; N];
        for elem in round.iter_mut() {
            *elem = read_field4(r)?;
        }
        out.push(round);
    }
    Ok(out)
}

fn write_vec_field(w: &mut impl Write, data: &[[u32; 4]]) -> io::Result<()> {
    for elem in data {
        write_u32_slice(w, elem)?;
    }
    Ok(())
}

fn read_vec_field(r: &mut impl Read, count: usize) -> io::Result<Vec<[u32; 4]>> {
    let mut out = Vec::with_capacity(count);
    for _ in 0..count {
        out.push(read_field4(r)?);
    }
    Ok(out)
}

fn read_field4(r: &mut impl Read) -> io::Result<[u32; 4]> {
    let mut buf = [0u32; 4];
    r.read_exact(bytemuck::cast_slice_mut::<u32, u8>(&mut buf))?;
    Ok(buf)
}
