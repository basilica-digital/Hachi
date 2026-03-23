mod fields;

use tfhe_ntt::prime32::Plan;

use std::slice;
use std::fs::File;
use std::io::{Write, BufWriter};
use std::alloc::{alloc, dealloc, Layout};
use std::ptr::NonNull;
use std::time::Instant;
use std::arch::x86_64::*;
use std::hint::black_box;

use rand::Rng;
use rand::distributions::{Distribution, Uniform};
use rand::thread_rng;
use rand::RngCore;

use rayon::prelude::*;

use crate::fields::fields::*;
include!("macro.rs");
include!("sumcheck.rs");
include!("rns.rs");
include!("foldwitness.rs");
include!("foldwitness_jolt.rs");
include!("commit.rs");


#[repr(C, align(64))]
#[derive(Clone, Copy)]
pub struct Align64<T>(pub T);

#[repr(align(64))]
#[derive(Clone)]
pub struct Align64U32(pub [u32; 16]); // 16 * 4 bytes = 64 bytes

pub unsafe fn ringmul(a: Vec<u32>, b: Vec<u32>) -> [u32; 1024] {
    let mut tmp = [0u64; 2048];
    let mut c = [0u32; 1024];
    
    let m = (1u64 << 32) - 99;

    for i in 0..1024 {
        for j in 0..1024 {
            let prod = a[i] as u64 * b[j] as u64;
            tmp[i + j] = (tmp[i + j] + prod) % m;
        }
    }

    for i in 0..1024 {
        if tmp[i] >= tmp[i + 1024] {
            c[i] = ((tmp[i] - tmp[i + 1024]) % m) as u32;
        } else {
            c[i] = ((tmp[i] + m - tmp[i + 1024]) % m) as u32;
        }
    }
    c
}

pub unsafe fn ringmul_jolt(z: &mut [i32], a: &[i32], b: &[u8]) {
    let mut tmp = [0i64; 2048];
    
    let m = (1i64 << 31) - 19;

    for i in 0..4 {
        let shift = b[i] as usize + i * 256;
        for j in 0..1024 {
            tmp[shift + j] += a[j] as i64;
        }
    }

    for i in 0..1024 {
        let res = (tmp[i] - tmp[i+1024]) % m;
        z[i] = if res < 0 { (res + m) as i32 } else { res as i32 };
    }
}

pub unsafe fn innerproduct_jolt(z: &mut [i32], c: &[i32], s: &[u8]){
    for x in z.iter_mut() { *x = 0; }
    let m = (1i64 << 31) - 19;
    
    for i in 0..(1<<16){
        let mut tmp_z = [0i32; 1024];
        ringmul_jolt(&mut tmp_z, &c[i*1024..(i+1)*1024], &s[i*4..(i+1)*4]);
        for j in 0..1024 {
            let res = (z[j] as i64 + tmp_z[j] as i64) % m;
            z[j] = if res < 0 { (res + m) as i32 } else { res as i32 };
        }
    }
}

pub fn ringmla_idx(z: &mut [i16], c: &[i16], s: &[i16]) {
    let mut tmp = [0i16; 2048];
    
    let m = (1i64 << 32) - 99;

    for i in 0..16 {
        if c[i] < 1024{
            for j in 0..1024{
                tmp[c[i] as usize + j] += s[j];
            }
        }else{
            for j in 0..1024{
                tmp[c[i] as usize - 1024 + j] -= s[j];
            }
        }
    }

    for i in 0..1024 {
        let res = tmp[i] - tmp[i+1024];
        z[i] += res;
    }
}


pub fn fold_witness_ref(z: &mut [i16], c: &[i16], s: &[i16]){
    for i in 0..(1<<13){
        for j in 0..(1<<10){
            ringmla_idx(& mut z[i<<10..(i+1)<<10], &c[j<<4..(j+1)<<4], &s[(i<<20)+(j<<10) ..(i<<20)+((j+1)<<10)]);
        }
    }
}


#[allow(non_snake_case)]
fn main()-> Result<(), Box<dyn std::error::Error>>{

    // -----> FOLDING WITNESS <-----
    let mut s: Vec<i16> = vec![0i16; 1<<33];
    thread_rng().fill(&mut s[..]);
    for val in s.iter_mut() { *val &= 0xF; }
    let c = generate_sparse_c_idx(1<<10);
    let mut z: Vec<i16> = vec![0; 1<<23];
    let mut zref: Vec<i16> = vec![0; 1<<23];

    let start = Instant::now();
    unsafe {
        fold_witness(&mut z, &c, &s);
    }
    let fw_duration = start.elapsed();
    println!("Folding witness: {:?}", fw_duration);
    // fold_witness_ref(&mut zref, &c, &s);
    // for i in 0..(1<<23){
    //     assert_eq!(z[i], zref[i]);
    // }
    // println!("Folding witness correct");


    // // -----> FOLDING WITNESS JOLT <-----
    // let height = 16;
    // let dimention = 10;
    // let width = 10;
    // let size_csv = format!("2^{}*2^{}", height, width);
    
    // // init s [height][width][dimention]
    // let n_s: usize = 1 << (height + width + 4);
    // let mut s = vec![0u8; n_s];
    // sparse_random_1_index(&mut s);

    // // init c [height][dimention]
    // let n_c: usize = 1 << (height + dimention);
    // let m = (1i64 << 31) - 19;
    // let mut c_buf = vec![Align64([0i32; 16]); n_c / 16];
    // let c = unsafe { std::slice::from_raw_parts_mut(c_buf.as_mut_ptr() as *mut i32, n_c) };
    // for val in c.iter_mut() {
    //     *val = thread_rng().gen_range(0..m as i32);
    // }

    // // init z [width][dimention]
    // let n_z: usize = 1 << (width + dimention);
    // let mut z_buf = vec![Align64([0i32; 16]); n_z / 16];
    // let z = unsafe { std::slice::from_raw_parts_mut(z_buf.as_mut_ptr() as *mut i32, n_z) };

    // // test z
    // let mut z_test = vec![0i32; n_z];

    // // folding witness
    // unsafe {
    //     let start = Instant::now();
    //     fold_witness_jolt_4(z, c, &s);
    //     let fwj_duration = start.elapsed();
    //     println!("CSV_RESULT: 4,{},{:?}", size_csv, fwj_duration);
    //     // fold_witness_jolt(z, c, &s, 1);
    //     // for i in 0..(n_z / 1024) {
    //     //     innerproduct_jolt(&mut z_test[i*1024..(i+1)*1024], &c, &s[i*((1<<18))..]);
    //     // }
    // }

    // // test z
    // for i in 0..n_z{
    //     assert_eq!(z[i], z_test[i]);
    // }


    // // -----> Commit <-----
    // let n = 1 << 10;
    // let height = 1 << 13;
    // let sz = 1 << 20;
    // let s = &mut generate_random_data_32bit(sz, n);
    // let acap = &mut generate_random_data_32bit(height, n);
    // let bcap = &mut generate_random_data_32bit(height, n);
    // let mut t = vec![Align64U32([0u32; 16]); 1<<19];
    // let mut u = Align64([0u32; 1024]);

    // let start = Instant::now();
    // unsafe{
    //     let t_ptr = t.as_mut_ptr() as *mut u32;
    //     let t_len = (1 << 19) * 16; 
    //     let t_flat_slice = std::slice::from_raw_parts_mut(t_ptr, t_len);
    //     black_box(commit(&mut u.0, t_flat_slice, s, acap, bcap));
    // }
    // let commit_duration = start.elapsed();
    // println!("Commit: {:?}", commit_duration);

    // // -----> Sumcheck <-----
    // sumcheck();


    Ok(())
}

#[cfg(test)]
mod ring_tests {
    use super::*;
    use std::arch::x86_64::*;
    unsafe fn ring_inner_product_ref(s: &[u32], a: &[u32], height: usize, n: usize) -> Vec<u32> {
        let mut acc = vec![0u64; n];
        let m = (1u64 << 32) - 99; // Goldilocks prime
        
        for j in 0..height {
            let s_poly = &s[j * n .. (j + 1) * n];
            let a_poly = &a[j * n .. (j + 1) * n];
            let prod = unsafe {ringmul(s_poly.to_vec(), a_poly.to_vec())}; 
            for i in 0..n {
                acc[i] = (acc[i] + prod[i] as u64) % m;
            }
        }
        acc.into_iter().map(|x| x as u32).collect()
    }

    #[test]
    fn test_ring_inner_product() {
        unsafe {
            let n = 1 << 10; // 1024
            let height = 1<<13; // 8192
            let moduli = [759207937, 759304193];
            
            let q1 = _mm512_set1_epi32(759207937);
            let q1_inv = _mm512_set1_epi32(754935809);
            let q2 = _mm512_set1_epi32(759304193);
            let q2_inv = _mm512_set1_epi32(331214849);
            
            let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();
            let random_ring = generate_random_data_32bit(height, n);
            let small_ring = generate_random_data_4bit(height, n);
            
            let mut u = vec![0u32; n];
            let u_ptr = u.as_mut_ptr() as *mut __m512i;
            
            let mut v1 = [_mm512_set1_epi32(0); 64];
            let mut v2 = [_mm512_set1_epi32(0); 64];

            for j in 0..height {
                let s_slice = &small_ring[j * n .. (j + 1) * n];
                let a_slice = &random_ring[j * n .. (j + 1) * n];

                let mut a1 = Align64([0u32; 1024]);
                let mut a2 = Align64([0u32; 1024]);
                let mut s1 = Align64([0u32; 1024]);
                let mut s2 = Align64([0u32; 1024]);

                // RNS Decompose A
                for k in (0..1024).step_by(16) {
                    rns_decompose_16_elements(
                        a_slice,
                        &mut a1.0,
                        &mut a2.0, 
                        k
                    );
                }
                s1.0.copy_from_slice(s_slice);
                s2.0.copy_from_slice(s_slice);
                plans[0].fwd(&mut a1.0);
                plans[1].fwd(&mut a2.0);
                plans[0].fwd(&mut s1.0);
                plans[1].fwd(&mut s2.0);

                let s1_ptr = s1.0.as_ptr() as *const __m512i;
                let s2_ptr = s2.0.as_ptr() as *const __m512i;
                let a1_ptr = a1.0.as_ptr() as *const __m512i;
                let a2_ptr = a2.0.as_ptr() as *const __m512i;

                for k in 0..64 {
                    let s1_vec = _mm512_load_si512(s1_ptr.add(k));
                    let a1_vec = _mm512_load_si512(a1_ptr.add(k));
                    let p1 = montproduct_759207937(s1_vec, a1_vec, q1, q1_inv);
                    
                    let s2_vec = _mm512_load_si512(s2_ptr.add(k));
                    let a2_vec = _mm512_load_si512(a2_ptr.add(k));
                    let p2 = montproduct_759304193(s2_vec, a2_vec, q2, q2_inv);
                    
                    let sum1 = _mm512_add_epi32(v1[k], p1);
                    let ge_q1 = _mm512_cmpge_epu32_mask(sum1, q1);
                    v1[k] = _mm512_mask_sub_epi32(sum1, ge_q1, sum1, q1);

                    let sum2 = _mm512_add_epi32(v2[k], p2);
                    let ge_q2 = _mm512_cmpge_epu32_mask(sum2, q2);
                    v2[k] = _mm512_mask_sub_epi32(sum2, ge_q2, sum2, q2);
                }
            }
            for k in 0..64 {
                v1[k] = barrett_mul_4194304_759207937(v1[k]);
                v2[k] = barrett_mul_4194304_759304193(v2[k]);
            }
            plans[0].inv(std::slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 1024));
            plans[1].inv(std::slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 1024));
            irns(&v1, &v2, u_ptr, 0);
            let m = (1u64 << 32) - 99;
            let u_ref = ring_inner_product_ref(&small_ring, &random_ring, height, n);

            let mut all_match = true;
            for i in 0..1024 {
                let expected = (u_ref[i] as u64) % m;
                let actual = (u[i] as u64) % m;
                if expected != actual {
                    println!("Failed at index {}: expected {}, got {}", i, expected, actual);
                    all_match = false;
                    break; 
                }
            }
            assert!(all_match, "Ring Inner Product test failed! The math does not match.");
            println!("Ring Inner Product test passed");
        }
    }
}