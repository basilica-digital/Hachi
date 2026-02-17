mod fields;

use tfhe_ntt::prime32::Plan;
use std::slice;
use std::fs::File;
use std::io::{Write, BufWriter};
use rand::Rng;
use std::time::Instant;
use std::arch::x86_64::*;
use std::hint::black_box;
use rand::distributions::{Distribution, Uniform};
use rand::thread_rng;
use rand::RngCore;

use rayon::prelude::*;

use crate::fields::fields::*;
include!("macro.rs");
include!("rns.rs");
include!("foldwitness.rs");
include!("commit.rs");

#[inline(never)]
fn wrapped_fwd(plan: &Plan, data: &mut [u32]) {
    plan.fwd(data);
}

#[inline(never)]
fn wrapped_inv(plan: &Plan, data: &mut [u32]) {
    plan.inv(data);
}

#[inline(always)]
pub unsafe fn add_ring(z: &mut [i16], s: &[i16]){
    for i in 0..1024 {
        z[i] = s[i];
    }
}

#[inline(always)]
pub unsafe fn addz(z: &mut [i16], s: &[i16]){
    for i in 0..1024{
        add_ring(z, &s[i*(1<<10)..]);
    }
}

#[inline(never)]
pub unsafe fn memory_test(z: &mut [i16], s: &[i16]) {
    z.par_chunks_mut(1 << 10).enumerate().for_each(|(i, z_row)| {
        addz(z_row, &s[i*(1 << 20)..]);
    });
}


#[inline(always)]
pub unsafe fn power(val: i16) -> i16{
    let mut v0 = val as i32;
    let mut v1 = val as i32 + 1;
    let mut v2 = val as i32 + 2;
    let mut v3 = val as i32 + 3;
    
    // 讓迴圈次數夠多，確保運算密集
    for _ in 0..200_000 {
        // 這四行指令互不相依，CPU 會嘗試在同一個週期內並行執行它們
        // 這會耗盡物理核心的所有 ALU 資源
        v0 = v0.wrapping_mul(v0).wrapping_add(12345); 
        v1 = v1.wrapping_mul(v1).wrapping_add(23456);
        v2 = v2.wrapping_mul(v2).wrapping_add(34567);
        v3 = v3.wrapping_mul(v3).wrapping_add(45678);
    }
    
    (v0 ^ v1 ^ v2 ^ v3) as i16
}



#[inline(never)]
pub unsafe fn core_test(a: &[i16]) -> [i16; 4096]{
    let mut b:[i16; 4096] = [0i16; 4096];
    b.par_iter_mut()
     .zip(a.par_iter()) 
     .for_each(|(dest, src)| {
         *dest = power(*src);
     });
    b
}

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



fn main()-> Result<(), Box<dyn std::error::Error>>{

    // // // -----> TEST NUM CORES <-----
    // let mut a:[i16; 4096] = [0i16; 4096];
    // for i in 0..4096{
    //     a[i] = i as i16+4096;
    // } 
    
    // let start = Instant::now();
    // unsafe{
    //     black_box(core_test(&a));
    // }
    
    // let duration = start.elapsed();
    // //println!("Folding witness: {:?}", fw_duration);
    // println!("Core test: {:?}", duration);


    // // -----> FOLDING WITNESS <-----
    // // init s [2^13][2^10][2^10]
    // let n: usize = 1 << 33;

    // let mut s: Vec<i16> = vec![0i16; n];
    // thread_rng().fill(&mut s[..]);
    // for val in s.iter_mut() { *val &= 0xF; }

    // // init c [2^10]
    // let n: usize = 1 << 20;
    // let mut c: Vec<i16> = (0..n).map(|_| thread_rng().gen_range(0..16)).collect();

    // // init z [2^13]
    // let n: usize = 1 << 23;
    // let mut z: Vec<i16> = vec![0; n];

    // let start = Instant::now();
    // // folding witness
    // unsafe {
    //     //fold_witness(&mut z, &c, &s);
    //     memory_test(&mut z, &s);
    // }
    // let fw_duration = start.elapsed();
    // //println!("Folding witness: {:?}", fw_duration);
    // println!("Memory test: {:?}", fw_duration);


    // -----> Commit <-----
    let start = Instant::now();
    let n = 1 << 10;
    let iterations = 1 << 10;
    let height = 1 << 13;
    let sz = 1 << 23;
    let s = &mut generate_random_data_4bit(sz, n);
    let acap1 = &mut generate_random_data_30bit(height, n);
    let acap2 = &mut generate_random_data_30bit(height, n);
    let mut u = vec![0u32; n];
    
    let init_duration = start.elapsed();

    println!("Init: {:?}", init_duration);
    let start = Instant::now();
    unsafe{
        commit(&mut u, s, acap1, acap2);
    }
    let commit_duration = start.elapsed();
    println!("Commit: {:?}", commit_duration);

    // // -----> Ring multiplication test <-----
    // // let file = File::create("out.txt")?;
    // // let mut writer = BufWriter::new(file);
    // unsafe{
    //     let n = 1 << 10;
    //     let moduli = [759207937, 759304193];
    //     let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();

    //     let random_ring = &mut generate_random_data_30bit(1, n);
    //     let small_ring = &mut generate_random_data_4bit(1, n);
        
    //     let mut r1:[u32;1024] = random_ring[0..1024].try_into().expect("Slice length mismatch");
    //     let mut r2:[u32;1024] = random_ring[0..1024].try_into().expect("Slice length mismatch");
    //     let mut s1:[u32;1024] = small_ring[0..1024].try_into().expect("Slice length mismatch");
    //     let mut s2:[u32;1024] = small_ring[0..1024].try_into().expect("Slice length mismatch");
    //     let mut u = vec![0u32; n];
    //     let u_ptr = u.as_mut_ptr() as *mut __m512i;
    //     let mut v1 = [_mm512_set1_epi32(0); 64];
    //     let mut v1_ptr = std::slice::from_raw_parts(v1.as_ptr() as *const u32, 1024);
    //     let mut v2 = [_mm512_set1_epi32(0); 64];
    //     let mut v2_ptr = std::slice::from_raw_parts(v2.as_ptr() as *const u32, 1024);
        
    //     plans[0].fwd(&mut s1);
    //     plans[1].fwd(&mut s2);
    //     plans[0].fwd(&mut r1);
    //     plans[1].fwd(&mut r2);
        
        
    //     let s1_ptr = s1.as_mut_ptr() as *mut __m512i;
    //     let s2_ptr = s2.as_mut_ptr() as *mut __m512i;
    //     let capa1_ptr = r1.as_mut_ptr() as *mut __m512i;
    //     let capa2_ptr = r2.as_mut_ptr() as *mut __m512i;
    //     for k in 0..64{
    //         let s1_vec = _mm512_loadu_si512(s1_ptr.add(k));
    //         let s2_vec = _mm512_loadu_si512(s2_ptr.add(k));
    //         let a1_vec = _mm512_loadu_si512(capa1_ptr.add(k));
    //         let a2_vec = _mm512_loadu_si512(capa2_ptr.add(k));
    //         v1[k] = simple_mul_mod_759207937(s1_vec, a1_vec);//montproduct_759207937(s1_vec, a1_vec);
    //         v2[k] = simple_mul_mod_759304193(s2_vec, a2_vec);//montproduct_759304193(s2_vec, a2_vec);
    //     }
        
    //     plans[0].inv(slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 1024));
    //     plans[1].inv(slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 1024));    
    //     irns(&v1, &v2, u_ptr, 0);

    //     let m = (1u64 << 32) - 99;

    //     let u_ref = ringmul(small_ring.to_vec(), random_ring.to_vec());

    //     for i in 0..1024 {
    //         if (u_ref[i] as u64) % m != (u[i] as u64) % m {
    //             println!("Failed at index {}: expected {}, got {}", i, (u_ref[i] as u64) % m, (u[i] as u64) % m);
    //         }
    //     }
    // }
    
    Ok(())
}