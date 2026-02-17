#[inline(always)]
pub unsafe fn rns_decompose_16_elements(
    input_data: &[u32],
    d1_target: &mut [u32],
    d2_target: &mut [u32],
    offset: usize,
) {
    let x = load_vec_to_m512(input_data, offset);
    let res1 = barrett_fake_759207937(x);
    let res2 = barrett_fake_759304193(x);
    store_m512_to_vec(res1, d1_target, offset);
    store_m512_to_vec(res2, d2_target, offset);
}


pub unsafe fn irns_vec(a: __m512i, b:__m512i) -> __m512i{
    let mut diff = _mm512_sub_epi32(b, a);
    let is_neg = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), diff);
    diff = _mm512_mask_add_epi32(diff, is_neg, diff, _mm512_set1_epi32(759304193));
	let v = barrett_mul_759304193(diff);
    mla_mod32(a, v)
}

pub unsafe fn irns(val1: &[__m512i; 64], val2: &[__m512i; 64], out: *mut __m512i, offset: usize){
    for i in 0..64{
        _mm512_storeu_si512(out.add(offset+i) as *mut _, irns_vec(val1[i], val2[i]));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::arch::x86_64::*;

    unsafe fn set_all_epi32(val: u32) -> __m512i {
        _mm512_set1_epi32(val as i32)
    }

    unsafe fn get_first_u32(vec: __m512i) -> u32 {
        let low_128 = _mm512_extracti32x4_epi32::<0>(vec);
        _mm_cvtsi128_si32(low_128) as u32
    }

    #[test]
    fn test_irns_vec_correctness() {
        if is_x86_feature_detected!("avx512f") {
            unsafe {
                let q1: u64 = 759207937;
                let q2: u64 = 759304193;
                
                
                // 測試測資 1
                let secret_val: u64 = 2885016803;
                let a = (secret_val % q1) as u32;
                let b = (secret_val % q2) as u32;

                let res_vec = irns_vec(set_all_epi32(a), set_all_epi32(b));
                let res = get_first_u32(res_vec) as u64;

                println!("q1:{}, q2:{}", res % q1, res % q2);

                // 驗證是否拼回原始值 (在 mod q1*q2 下)
                // 如果 secret_val < 2^32，res 應該直接等於 secret_val
                assert_eq!(res % q1, a as u64, "q1 剩餘不符");
                assert_eq!(res % q2, b as u64, "q2 剩餘不符");

                // 大規模隨機測試
                let mut seed: u32 = 0xABCDE;
                for _ in 0..100000 {
                    seed = seed.wrapping_mul(1103515245).wrapping_add(12345);
                    let test_val = seed as u64; // 隨機產生一個 32-bit 範圍內的數
                    
                    let a_simd = (test_val % q1) as u32;
                    let b_simd = (test_val % q2) as u32;
                    
                    let v_res = irns_vec(set_all_epi32(a_simd), set_all_epi32(b_simd));
                    let actual = get_first_u32(v_res) as u64;
                    
                    if actual % q1 != a_simd as u64 || actual % q2 != b_simd as u64 {
                        panic!("驗證失敗！原始值: {}, 得到的 a: {}, b: {}, 恢復值: {}", 
                               test_val, a_simd, b_simd, actual);
                    }
                }
                println!("✅ 100,000 組 IRNS 隨機測試通過！");
            }
        }
    }
}