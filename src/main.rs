mod fields;

use tfhe_ntt::prime32::Plan;
use std::slice;
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



fn main(){

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


    let n = 1024;
    // let moduli = [817153, 833537, 854017, 759207937, 759304193];
    let moduli = [759207937, 759304193];
    let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();

    // let mut data_sets: Vec<Vec<u32>> = generate_random_data(n);

    // let mut seed = 0x12345678u64;
    // let mut next_rand = || {
    //     seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    //     seed
    // };

    // for (i, &q) in moduli.iter().enumerate() {
    //     let mut all_tests_passed = true;

    //     for test_idx in 0..10000 {
    //         // 生成隨機測資，確保每個數都在 [0, q-1]
    //         let original_data: Vec<u32> = (0..n).map(|_| (next_rand() % q as u64) as u32).collect();
    //         let mut test_data = original_data.clone();

    //         // 正逆轉換
    //         plans[i].fwd(&mut test_data);
    //         plans[i].inv(&mut test_data);

    //         // 歸一化 (除以 N)
    //         let n_inv = mod_inverse(n as u64, q as u64) as u32;
    //         for val in test_data.iter_mut() {
    //             *val = ((*val as u64 * n_inv as u64) % q as u64) as u32;
    //         }

    //         // 嚴格比對
    //         if original_data != test_data {
    //             println!("❌ 模數 {} 驗證失敗於隨機組 #{}！", q, test_idx);
    //             all_tests_passed = false;
    //             break;
    //         }
    //     }

    //     if all_tests_passed {
    //         println!("✅ 模數 {} (約 {} bit) 通過 1000 組隨機壓力測試", q, 32 - q.leading_zeros());
    //     }
    // }

    // -----> RNS + NTT test <-----

    // --- Warm-up ---
    // let sz = 100;
    // let mut data_sets: Vec<u32> = generate_random_data(sz, n);
    // for zz in 0..100 {
    //     plans[0].fwd(&mut data_sets[zz*1024..(zz+1)*1024]);
    // }

    // --- (20-bit x3) ---
    // let start_small = Instant::now();
    // for _ in 0..iterations {
    //     for i in 0..3 {
    //         plans[i].fwd(&mut data_sets[0]);
    //     }
    // }
    // let duration_small = start_small.elapsed();

    // --- (30-bit x2) ---
    println!("Start init");
    let start = Instant::now();
    let n = 1 << 10;
    let iterations = 1 << 10;
    let height = 1 << 13;
    let sz = 1 << 23;

    let s = &mut generate_random_data_4bit(sz, n);
    let acap = &mut generate_random_data_30bit(height, n);
    let bcap = &mut generate_random_data_30bit(height, n);
    let mut t = vec![0u32; iterations*n];
    let mut u = vec![0u32; n];
    let mut gt = vec![0u32; height*n];
    let mut capa1 = vec![0u32; height*n];
    let mut capa2 = vec![0u32; height*n];
    let mut capb1 = vec![0u32; height*n];
    let mut capb2 = vec![0u32; height*n];


    let gt_ptr = gt.as_mut_ptr() as *mut __m512i;
    let t_ptr = t.as_mut_ptr() as *mut __m512i;
    let u_ptr = u.as_mut_ptr() as *mut __m512i;
    let capa1_ptr = capa1.as_mut_ptr() as *mut __m512i;
    let capa2_ptr = capa2.as_mut_ptr() as *mut __m512i;
    let capb1_ptr = capb1.as_mut_ptr() as *mut __m512i;
    let capb2_ptr = capb2.as_mut_ptr() as *mut __m512i;
    let init_duration = start.elapsed();
    println!("Init: {:?}", init_duration);
    
    
    let start = Instant::now();
    unsafe{
        // RNS on A
        for i in 0..height{
            for j in (0..1024).step_by(16) {
                rns_decompose_16_elements(&acap[i*1024+j..i*1024+(j+16)], &mut capa1[i*1024+j..i*1024+(j+16)], &mut capa2[i*1024+j..i*1024+(j+16)], j);
            }
        }

        // As
        for i in 0..iterations{
            let mut val1 = [_mm512_set1_epi32(0); 64];
            let mut val2 = [_mm512_set1_epi32(0); 64];
            for j in 0..height{
                let mut s1:[u32;1024] = (&s[j*1024+i..(j+1)*1024+i]).try_into().expect("Slice length mismatch");
                let mut s2:[u32;1024] = (&s[j*1024+i..(j+1)*1024+i]).try_into().expect("Slice length mismatch");
                wrapped_fwd(&plans[0], &mut s1);
                wrapped_fwd(&plans[1], &mut s2);
                let s1_ptr = s1.as_mut_ptr() as *mut __m512i;
                let s2_ptr = s2.as_mut_ptr() as *mut __m512i;
                for k in 0..64{
                    let s1_vec = _mm512_loadu_si512(s1_ptr.add(k));
                    let s2_vec = _mm512_loadu_si512(s2_ptr.add(k));
                    let a1_vec = _mm512_loadu_si512(capa1_ptr.add(k));
                    let a2_vec = _mm512_loadu_si512(capa2_ptr.add(k));
                    val1[k] = barrett_fake_759207937(_mm512_add_epi32(val1[k], montproduct_759207937(s1_vec, a1_vec)));
                    val2[k] = barrett_fake_759304193(_mm512_add_epi32(val2[k], montproduct_759304193(s2_vec, a2_vec)));  
                }
            }
            wrapped_inv(&plans[0], slice::from_raw_parts_mut(val1.as_mut_ptr() as *mut u32, 1024));
            wrapped_inv(&plans[1], slice::from_raw_parts_mut(val2.as_mut_ptr() as *mut u32, 1024));
            irns(&val1, &val2, t_ptr, i);
        }

        // decompose t
        let mask15 = _mm512_set1_epi32(15);
        for i in 0..iterations{
            for j in 0..64{
                let mut val = _mm512_loadu_si512(t_ptr.add(i*64+j));
                for k in 0..8{
                    _mm512_storeu_si512((gt_ptr.add(i*512+k*64+j)) as *mut _ , _mm512_and_epi32(mask15, val));
                    val = _mm512_srli_epi32(val, 4);
                }
            }
        }
        // RNS on B
        for i in 0..height{
            for j in (0..1024).step_by(16) {
                rns_decompose_16_elements(&bcap[i*1024+j..i*1024+(j+16)], &mut capb1[i*1024+j..i*1024+(j+16)], &mut capb2[i*1024+j..i*1024+(j+16)], j);
            }
        }
        // B*gt
        let mut v1 = [_mm512_set1_epi32(0); 64];
        let mut v2 = [_mm512_set1_epi32(0); 64];
        for i in 0..height{
            let mut gt1:[u32;1024] = (&gt[i*1024..(i+1)*1024]).try_into().expect("Slice length mismatch");
            let mut gt2:[u32;1024] = (&gt[i*1024..(i+1)*1024]).try_into().expect("Slice length mismatch");
            wrapped_fwd(&plans[0], &mut gt1);
            wrapped_fwd(&plans[1], &mut gt2);
            let gt1_ptr = gt1.as_mut_ptr() as *mut __m512i;
            let gt2_ptr = gt2.as_mut_ptr() as *mut __m512i;
            for j in 0..64{
                let gt1_vec = _mm512_loadu_si512(gt1_ptr.add(j));
                let gt2_vec = _mm512_loadu_si512(gt2_ptr.add(j));
                let b1_vec = _mm512_loadu_si512(capb1_ptr.add(j));
                let b2_vec = _mm512_loadu_si512(capb2_ptr.add(j));
                v1[j] = barrett_fake_759207937(_mm512_add_epi32(v1[j], montproduct_759207937(gt1_vec, b1_vec)));
                v2[j] = barrett_fake_759304193(_mm512_add_epi32(v2[j], montproduct_759304193(gt2_vec, b2_vec)));  
            }
        }
        wrapped_inv(&plans[0], slice::from_raw_parts_mut(v1.as_mut_ptr() as *mut u32, 1024));
        wrapped_inv(&plans[1], slice::from_raw_parts_mut(v2.as_mut_ptr() as *mut u32, 1024));
        irns(&v1, &v2, u_ptr, 0);
        black_box(u);
    }
    let commit_duration = start.elapsed();
    println!("Commit: {:?}", commit_duration);

}









    
        
        
    //     // unsafe {
    //     //     for i in (0..n).step_by(16) {
    //     //         rns_decompose_16_elements(current_data, &mut d1, &mut d2, i);
    //     //     }
    //     // }
    //     rns_duration += start.elapsed();
    //     let start = Instant::now();
    //     wrapped_fwd(&plans[0], &mut d1);
    //     wrapped_fwd(&plans[1], &mut d2);
    //     ntt_duration += start.elapsed();
    //     black_box(&d1);
    //     black_box(&d2);
    // }

    //println!("Repeated {} iterations", iterations);
    // println!("前三組平均單次 NTT 耗時: {:?}", duration_small / (iterations * 3));
    // println!("RNS: {:?}", rns_duration);
    // println!("NTT: {:?}", ntt_duration);
    // println!("後兩組平均單次 NTT 耗時: {:?}", total_duration / (iterations * 2));




//     // // 3. Dot Product (Hadamard Product)
//     // let start = Instant::now();
//     // for i in 0..n {
//     //     a[i] = ((a[i] as u64 * b[i] as u64) % q as u64) as u32;
//     // }

//     // 4. Backward NTT (使用 inv)
//     let start = Instant::now();
//     plan.inv(&mut a);

//     let scale = n as u32;
//     println!("驗證第一項: 預期 {}, 實際 {}", (data_backup[0] * scale) % q, a[0] % q);
// }
