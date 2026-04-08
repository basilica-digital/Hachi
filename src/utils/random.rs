use crate::{Align64, AlignedU32Vec};
use rand::prelude::*;
use rayon::prelude::*;
use rand::distributions::{Distribution, Uniform};
use crate::utils::ds::*;
use rand_chacha::ChaCha12Rng;

pub fn generate_random_data_4bit_packed(sz: usize, n: usize) -> AlignedU8Vec {
    let total_elements = sz * n;
    assert!(total_elements % 32 == 0, "Total length must be a multiple of 32"); 
    let total_bytes = total_elements / 2;
    let mut data = AlignedU8Vec {
        inner: vec![Align64([0u8; 64]); total_bytes / 64 + 1],
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

// pub fn generate_sparse_c_idx(rng: &mut ChaCha12Rng, num_rings: usize) -> Vec<i16> {
//     let n = 16;
//     let mut c = Vec::with_capacity(num_rings * n);

//     for _ in 0..num_rings{
//         let mut group: Vec<i16> = rand::seq::index::sample(rng, 1024, n)
//             .into_iter()
//             .map(|x| x as i16)
//             .collect();
            
//         group.sort_unstable();
        
//         let num_to_offset = rng.gen_range(0..=n);
//         let mut indices: Vec<usize> = (0..n).collect();
//         indices.shuffle(rng);
        
//         for &idx in &indices[..num_to_offset] {
//             group[idx] += 1024;
//         }
//         c.extend(group);
//     }
//     c
// }