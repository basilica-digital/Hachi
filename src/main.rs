#![allow(warnings)]
pub mod field;
pub mod ntt;

use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::ntt::ntt::*;
use crate::ntt::intt::*;
use crate::ntt::transpose::*;

fn main() {
    println!("Hello, world!");
}

// #[cfg(test)]
mod tests {
    use std::arch::x86_64::*;
    use crate::field::fields::*;
    use crate::ntt::ntt::*;
    use crate::ntt::intt::*;
    use crate::ntt::transpose::*;

    unsafe fn load_i16_to_m256(input: &[i16; 256]) -> [__m256i; 16] {
        let mut output: [__m256i; 16] = [_mm256_setzero_si256(); 16];
        for i in 0..16 {
            output[i] = _mm256_loadu_si256(input.as_ptr().add(i * 16) as *const __m256i);
        }
        output
    }

    macro_rules! test_ntt_roundtrip {
        ($func_name:ident, $q:expr, $ntt_fn:ident, $intt_fn:ident) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }
                
                unsafe {
                    let q = $q as i32;
                    let scaling_factor = 256i32;

                    let mut p_data = [0i16; 256];
                    for i in 0..256 {
                        p_data[i] = (700*i as i32 % q) as i16;
                    }
                    let p_original = p_data.clone();

                    let ptr = p_data.as_ptr() as *const [__m256i; 16];
                    let input_vecs = std::ptr::read_unaligned(ptr);

                    let ntt_res = crate::ntt::ntt::$ntt_fn(input_vecs);
                    let intt_res = crate::ntt::intt::$intt_fn(ntt_res);

                    let mut p_actual = [0i16; 256];
                    let out_ptr = p_actual.as_mut_ptr() as *mut [__m256i; 16];
                    std::ptr::write_unaligned(out_ptr, intt_res);

                    for i in 0..256 {
                        let a = p_original[i] as i32;
                        let b = p_actual[i] as i32;
                        
                        let expected_mod_q = (a * scaling_factor).rem_euclid(q);
                        let actual_mod_q = b.rem_euclid(q);

                        let mut error_indices = Vec::new();
                        for i in 0..256 {
                            let a = p_original[i] as i32;
                            let b = p_actual[i] as i32;
                            
                            let expected_mod_q = (a * scaling_factor).rem_euclid(q);
                            let actual_mod_q = b.rem_euclid(q);

                            if actual_mod_q != expected_mod_q {
                                error_indices.push(i);
                            }
                        }
                        if !error_indices.is_empty() {
                            println!("\n=== DEBUG: Round-trip FAILED (Congruence Check) for q={} ===", q);
                            println!("Total errors found: {}/256", error_indices.len());
                            println!("{:<10} | {:<10} | {:<10} | {:<10}", "Index", "Original", "Expected", "Actual(Raw)");
                            println!("{:-<50}", "");

                            for &i in &error_indices {
                                let a = p_original[i] as i32;
                                let mut b = p_actual[i] as i32;
                                let expected = (a * scaling_factor).rem_euclid(q);
                                
                                while b < 0{
                                    b += q;
                                }
                                println!("{:<10} | {:<10} | {:<10} | {:<10}", i, a, expected, b);
                            }
                            println!("==========================================\n");

                            panic!(
                                "Round-trip FAILED for q={}: {} mismatches found.", 
                                q, error_indices.len()
                            );
                        }
                    }
                    println!("SUCCESS: Round-trip for q={} passed (Congruent within modulo).", q);
                }
            }
        };
    }
    test_ntt_roundtrip!(roundtrip_7681, 7681, ntt_7681, intt_7681);
    test_ntt_roundtrip!(roundtrip_10753, 10753, ntt_10753, intt_10753);
    test_ntt_roundtrip!(roundtrip_11777, 11777, ntt_11777, intt_11777);
    test_ntt_roundtrip!(roundtrip_12289, 12289, ntt_12289, intt_12289);
    test_ntt_roundtrip!(roundtrip_13313, 13313, ntt_13313, intt_13313);
    test_ntt_roundtrip!(roundtrip_15361, 15361, ntt_15361, intt_15361);
}
