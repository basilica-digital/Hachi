#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

use std::arch::x86_64::*;

#[target_feature(enable = "avx512f")]
pub unsafe fn _mm512_mulhi_epu32(a: __m512i, b: __m512i) -> __m512i {
    let mul_evens = _mm512_mul_epu32(a, b);
    let a_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(a);
    let b_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(b);
    let mul_odds = _mm512_mul_epu32(a_odds, b_odds);
    let hi_evens = _mm512_srli_epi64::<32>(mul_evens);
    let mask = 0xAAAA; 
    _mm512_mask_blend_epi32(mask, hi_evens, mul_odds)
}

// barrett_fake
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn barrett_fake_759207937(x: __m512i) -> __m512i{
    let d = _mm512_mulhi_epu32(x, _mm512_set1_epi32(6));
    let e = _mm512_mullo_epi32(d, _mm512_set1_epi32(759207937));
    _mm512_sub_epi32(_mm512_add_epi32(x, _mm512_set1_epi32(759207937)), e)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn barrett_fake_759304193(x: __m512i) -> __m512i{
    let d = _mm512_mulhi_epu32(x, _mm512_set1_epi32(6));
    let e = _mm512_mullo_epi32(d, _mm512_set1_epi32(759304193));
    _mm512_sub_epi32(_mm512_add_epi32(x, _mm512_set1_epi32(759304193)), e)
}

#[cfg(test)]
mod tests {
    use super::*;
    use core::arch::x86_64::*;

    unsafe fn dump_m512(v: __m512i) -> [u32; 16] {
        let mut arr = [0u32; 16];
        _mm512_storeu_si512(arr.as_mut_ptr() as *mut u32 as *mut __m512i, v);
        arr
    }

    macro_rules! test_barrett {
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx512f") { return; }
                use rand::Rng;

                unsafe {
                    let q = $q as u64;
                    let mut rng = rand::thread_rng();
                    let mut inputs = [0u32; 16];
                    
                    for _ in 0..1000000 {
                        for i in 0..16 {
                            inputs[i] = rng.gen::<u32>();
                        }

                        let input_vec = _mm512_loadu_si512(inputs.as_ptr() as *const __m512i);
                        let output_vec = super::$func_name(input_vec);
                        let outputs = dump_m512(output_vec);

                        for i in 0..16 {
                            assert_eq!(
                                inputs[i] as u64 % q, 
                                outputs[i] as u64 % q, 
                                "Fail at index {}, input: {}, q: {}", i, inputs[i], q
                            );
                        }
                    }
                }
            }
        };
    }

    test_barrett!(barrett_fake_759207937, 759207937);
    test_barrett!(barrett_fake_759304193, 759304193);
}