use crate::{Align64, AlignedU32Vec};
use rand::prelude::*;
use rayon::prelude::*;
use rand::distributions::{Distribution, Uniform};
use crate::utils::ds::*;
use rand_chacha::ChaCha12Rng;
use rand::SeedableRng;
use std::fs::File;
use std::io::{self, BufWriter, Write};
use std::path::Path;

pub fn generate_random_data_4bit_packed(sz: usize, n: usize) -> AlignedU8Vec {
    let total_elements = sz * n;
    assert!(total_elements % 32 == 0, "Total length must be a multiple of 32");
    let total_bytes = total_elements / 2;
    generate_random_bytes_4bit_packed(total_bytes)
}

/// Stream `total_bytes` of random 4-bit-packed data directly to `path`.
///
/// Memory usage is bounded to `chunk_bytes` regardless of `total_bytes`, so
/// this is the path for witnesses larger than RAM. Each chunk is filled in
/// parallel (Rayon) and then written sequentially via a `BufWriter`.
///
/// Both `total_bytes` and `chunk_bytes` must be multiples of 16.
pub fn stream_random_witness_to(
    path: &Path,
    total_bytes: u64,
    chunk_bytes: usize,
) -> io::Result<()> {
    assert!(chunk_bytes > 0 && chunk_bytes % 16 == 0);
    assert!(total_bytes % 16 == 0);

    let f = File::create(path)?;
    // 8 MiB kernel-side write buffer: big enough to keep NVMe fed without
    // adding meaningful RAM overhead on top of `chunk_bytes`.
    let mut w = BufWriter::with_capacity(8 * 1024 * 1024, f);

    let mut buf = vec![0u8; chunk_bytes];
    let mut remaining = total_bytes;
    while remaining > 0 {
        let this_chunk = remaining.min(chunk_bytes as u64) as usize;
        let slice = &mut buf[..this_chunk];
        slice.par_chunks_exact_mut(16).for_each(|c| {
            let mut rng = rand::thread_rng();
            for i in 0..16 {
                let low = (rng.next_u32() & 0xF) as u8;
                let high = (rng.next_u32() & 0xF) as u8;
                c[i] = (high << 4) | low;
            }
        });
        w.write_all(slice)?;
        remaining -= this_chunk as u64;
    }
    w.flush()?;
    Ok(())
}

/// Generate `total_bytes` of random 4-bit-packed data (two nibbles per byte).
///
/// `total_bytes` must be a multiple of 16 — the SIMD-friendly chunk size used
/// by the parallel fill loop. The backing store is 64-byte aligned and may
/// contain up to 63 bytes of zero padding past `len` for kernel safety.
pub fn generate_random_bytes_4bit_packed(total_bytes: usize) -> AlignedU8Vec {
    assert!(
        total_bytes % 16 == 0,
        "total_bytes must be a multiple of 16 (got {total_bytes})"
    );
    let num_chunks = total_bytes.div_ceil(64).max(1);
    let mut data = AlignedU8Vec {
        inner: vec![Align64([0u8; 64]); num_chunks],
        len: total_bytes,
    };

    data.par_chunks_exact_mut(16).for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for i in 0..16 {
            let low = (rng.next_u32() & 0xF) as u8;
            let high = (rng.next_u32() & 0xF) as u8;
            chunk[i] = (high << 4) | low;
        }
    });
    data
}

pub fn generate_random_q_element(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    data.par_iter_mut().for_each(|chunk| {
        let mut rng = rand::thread_rng();
        for x in chunk.0.iter_mut() {
            *x = rng.gen_range(0..(u32::MAX - 98));
        }
    });
    AlignedU32Vec { inner: data, len: total_u32 }
}

/// Deterministic variant of [`generate_random_q_element`].
///
/// Draws a per-call seed from `master_rng`, then derives independent
/// per-chunk sub-RNGs via `blake3::derive_key` so that Rayon's parallel
/// iteration order does not affect the output.
pub fn generate_random_q_element_seeded(
    sz: usize,
    n: usize,
    master_rng: &mut ChaCha12Rng,
) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;

    let mut call_seed = [0u8; 32];
    master_rng.fill_bytes(&mut call_seed);

    let mut data = vec![Align64([0u32; 16]); total_chunks];
    data.par_iter_mut().enumerate().for_each(|(i, chunk)| {
        let mut key_material = [0u8; 40]; // 32 (call_seed) + 8 (index)
        key_material[..32].copy_from_slice(&call_seed);
        key_material[32..].copy_from_slice(&i.to_le_bytes());
        let chunk_seed = blake3::derive_key("hachi-setup-chunk", &key_material);
        let mut rng = ChaCha12Rng::from_seed(chunk_seed);
        for x in chunk.0.iter_mut() {
            *x = rng.gen_range(0..(u32::MAX - 98));
        }
    });
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_fs_eval_points_q(rng: &mut ChaCha12Rng, sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    let ring_size = n / 16;
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    let q = 4294967197u32;
    for i in 0..sz {
        data[i * ring_size].0[0] = rng.gen_range(0..q);
    }
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_random_eval_points_q(sz: usize, n: usize) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    let ring_size = n / 16;
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    let q = 4294967197u32;
    let mut rng = rand::thread_rng();
    for i in 0..sz {
        data[i * ring_size].0[0] = rng.gen_range(0..q);
    }
    AlignedU32Vec { inner: data, len: total_u32 }
}

/// Deterministic variant of [`generate_random_eval_points_q`].
pub fn generate_random_eval_points_q_seeded(
    sz: usize,
    n: usize,
    rng: &mut ChaCha12Rng,
) -> AlignedU32Vec {
    let total_u32 = sz * n;
    assert!(total_u32 % 16 == 0, "Total length must be a multiple of 16");
    let total_chunks = total_u32 / 16;
    let ring_size = n / 16;
    let mut data = vec![Align64([0u32; 16]); total_chunks];
    let q = 4294967197u32;
    for i in 0..sz {
        data[i * ring_size].0[0] = rng.gen_range(0..q);
    }
    AlignedU32Vec { inner: data, len: total_u32 }
}

pub fn generate_sparse_c_idx(num_rings: usize) -> Vec<i16> {
    let n = 16;
    let mut c = Vec::with_capacity(num_rings * n);
    let mut rng = rand::thread_rng();
    for _ in 0..num_rings {
        let mut group: Vec<i16> = rand::seq::index::sample(&mut rng, 1024, n)
            .into_iter()
            .map(|x| x as i16)
            .collect();
        group.sort_unstable();
        let num_to_offset = rng.gen_range(0..=n);
        let mut indices: Vec<usize> = (0..n).collect();
        indices.shuffle(&mut rng);
        for &idx in &indices[..num_to_offset] {
            group[idx] += 1024;
        }
        c.extend(group);
    }
    c
}

/// Deterministic variant of [`generate_sparse_c_idx`].
pub fn generate_sparse_c_idx_seeded(num_rings: usize, rng: &mut ChaCha12Rng) -> Vec<i16> {
    let n = 16;
    let mut c = Vec::with_capacity(num_rings * n);
    for _ in 0..num_rings {
        let mut group: Vec<i16> = rand::seq::index::sample(rng, 1024, n)
            .into_iter()
            .map(|x| x as i16)
            .collect();
        group.sort_unstable();
        let num_to_offset = rng.gen_range(0..=n);
        let mut indices: Vec<usize> = (0..n).collect();
        indices.shuffle(rng);
        for &idx in &indices[..num_to_offset] {
            group[idx] += 1024;
        }
        c.extend(group);
    }
    c
}