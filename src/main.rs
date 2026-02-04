mod fields;

use tfhe_ntt::prime32::Plan;
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


fn main(){

    // -----> FOLDING WITNESS <-----
    // init s [2^13][2^10][2^10]
    let n: usize = 1 << 33;

    let mut s: Vec<i16> = vec![0i16; n];
    thread_rng().fill(&mut s[..]);
    for val in s.iter_mut() { *val &= 0xF; }

    // init c [2^10]
    let n: usize = 1 << 20;
    let mut c: Vec<i16> = (0..n).map(|_| thread_rng().gen_range(0..16)).collect();

    // init z [2^13]
    let n: usize = 1 << 23;
    let mut z: Vec<i16> = vec![0; n];

    let start = Instant::now();
    // folding witness
    unsafe {
        fold_witness(&mut z, &c, &s);
    }
    let fw_duration = start.elapsed();
    println!("Folding witness: {:?}", fw_duration);


    // let n = 1024;
    // // let moduli = [817153, 833537, 854017, 759207937, 759304193];
    // let moduli = [759207937, 759304193];
    // let plans: Vec<Plan> = moduli.iter().map(|&q| Plan::try_new(n, q).unwrap()).collect();

    // // let mut data_sets: Vec<Vec<u32>> = generate_random_data(n);

    // // let mut seed = 0x12345678u64;
    // // let mut next_rand = || {
    // //     seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
    // //     seed
    // // };

    // // for (i, &q) in moduli.iter().enumerate() {
    // //     let mut all_tests_passed = true;

    // //     for test_idx in 0..10000 {
    // //         // 生成隨機測資，確保每個數都在 [0, q-1]
    // //         let original_data: Vec<u32> = (0..n).map(|_| (next_rand() % q as u64) as u32).collect();
    // //         let mut test_data = original_data.clone();

    // //         // 正逆轉換
    // //         plans[i].fwd(&mut test_data);
    // //         plans[i].inv(&mut test_data);

    // //         // 歸一化 (除以 N)
    // //         let n_inv = mod_inverse(n as u64, q as u64) as u32;
    // //         for val in test_data.iter_mut() {
    // //             *val = ((*val as u64 * n_inv as u64) % q as u64) as u32;
    // //         }

    // //         // 嚴格比對
    // //         if original_data != test_data {
    // //             println!("❌ 模數 {} 驗證失敗於隨機組 #{}！", q, test_idx);
    // //             all_tests_passed = false;
    // //             break;
    // //         }
    // //     }

    // //     if all_tests_passed {
    // //         println!("✅ 模數 {} (約 {} bit) 通過 1000 組隨機壓力測試", q, 32 - q.leading_zeros());
    // //     }
    // // }

    // // -----> RNS + NTT test <-----

    // // --- Warm-up ---
    // // let sz = 100;
    // // let mut data_sets: Vec<u32> = generate_random_data(sz, n);
    // // for zz in 0..100 {
    // //     plans[0].fwd(&mut data_sets[zz*1024..(zz+1)*1024]);
    // // }

    // // --- (20-bit x3) ---
    // // let start_small = Instant::now();
    // // for _ in 0..iterations {
    // //     for i in 0..3 {
    // //         plans[i].fwd(&mut data_sets[0]);
    // //     }
    // // }
    // // let duration_small = start_small.elapsed();

    // // --- (30-bit x2) ---
    // let iterations = 1 << 23;
    // let mut d1 = vec![0u32; n];
    // let mut d2 = vec![0u32; n];

    // let mut rns_duration = std::time::Duration::new(0, 0);
    // let mut ntt_duration = std::time::Duration::new(0, 0);
    
    // let sz = 1 << 23;
    // let original_data = &generate_random_data(sz, n);
    // for zz in 0..iterations {
    //     let current_data = &original_data[zz*1024..(zz+1)*1024];
    //     let start = Instant::now();
    //     unsafe {
    //         for i in (0..n).step_by(16) {
    //             rns_decompose_16_elements(current_data, &mut d1, &mut d2, i);
    //         }
    //     }
    //     rns_duration += start.elapsed();
    //     let start = Instant::now();
    //     wrapped_fwd(&plans[0], &mut d1);
    //     wrapped_fwd(&plans[1], &mut d2);
    //     ntt_duration += start.elapsed();
    //     black_box(&d1);
    //     black_box(&d2);
    // }

    // println!("Repeated {} iterations", iterations);
    // // println!("前三組平均單次 NTT 耗時: {:?}", duration_small / (iterations * 3));
    // println!("RNS: {:?}", rns_duration);
    // println!("NTT: {:?}", ntt_duration);
    // // println!("後兩組平均單次 NTT 耗時: {:?}", total_duration / (iterations * 2));
}



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
