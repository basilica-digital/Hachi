#[cfg(target_arch = "x86_64")]
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

#[target_feature(enable = "avx512f")]
pub unsafe fn _mm512_mulhi_epu32_const(a: __m512i, b: __m512i) -> __m512i {
    let mul_evens = _mm512_mul_epu32(a, b);
    let a_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(a);
    let mul_odds = _mm512_mul_epu32(a_odds, b);
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

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
// montgomery multiplication
pub unsafe fn montproduct_759207937(x: __m512i, y: __m512i, q: __m512i, q_inv: __m512i) -> __m512i{
    let lo = _mm512_mullo_epi32(x, y);
    let hi = _mm512_mulhi_epu32(x, y);
    let d = _mm512_mullo_epi32(lo, q_inv);
    let e = _mm512_mulhi_epu32_const(d, q);
    let res = _mm512_sub_epi32(hi, e);
    let is_negative = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), res);
    _mm512_mask_add_epi32(res, is_negative, res, q)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn montproduct_759304193(x: __m512i, y: __m512i, q: __m512i, q_inv: __m512i) -> __m512i{
    let lo = _mm512_mullo_epi32(x, y);
    let hi = _mm512_mulhi_epu32(x, y);
    let d = _mm512_mullo_epi32(lo, q_inv);
    let e = _mm512_mulhi_epu32_const(d, q);
   let res = _mm512_sub_epi32(hi, e);
    let is_negative = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), res);
    _mm512_mask_add_epi32(res, is_negative, res, q)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn simple_mul_mod_759207937(x: __m512i, y: __m512i) -> __m512i {
    let q = 759207937i64;
    let mut x_arr = [0u32; 16];
    let mut y_arr = [0u32; 16];
    let mut res_arr = [0u32; 16];
    
    _mm512_storeu_si512(x_arr.as_mut_ptr() as *mut _, x);
    _mm512_storeu_si512(y_arr.as_mut_ptr() as *mut _, y);

    for i in 0..16 {
        let mut product = (x_arr[i] as u64) * (y_arr[i] as u64);
        product = product % q as u64;
        product = (product as u64) * (758466523 as u64);
        res_arr[i] = (product % q as u64) as u32;
    }

    _mm512_loadu_si512(res_arr.as_ptr() as *const _)
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn simple_mul_mod_759304193(x: __m512i, y: __m512i) -> __m512i {
    let q = 759304193i64;
    let mut x_arr = [0u32; 16];
    let mut y_arr = [0u32; 16];
    let mut res_arr = [0u32; 16];
    
    _mm512_storeu_si512(x_arr.as_mut_ptr() as *mut _, x);
    _mm512_storeu_si512(y_arr.as_mut_ptr() as *mut _, y);

    for i in 0..16 {
        let mut product = (x_arr[i] as u64) * (y_arr[i] as u64);
        product = product % q as u64;
        product = (product as u64) * (758562685 as u64);
        res_arr[i] = (product % q as u64) as u32;
    }

    _mm512_loadu_si512(res_arr.as_ptr() as *const _)
}

// a*161561972 % 759304193
#[inline(always)]
pub unsafe fn barrett_mul_759304193(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(759304193);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(913867449));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(161561972));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-759304193));
    let v = _mm512_add_epi32(d, g);
    let exceeds_q = _mm512_cmpge_epu32_mask(v, q);
    _mm512_mask_sub_epi32(v, exceeds_q, v, q)
}

#[inline(always)]
pub unsafe fn barrett_mul_4194304_759207937(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(759207937);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(23727884));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(4194304));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-759207937));
    let v = _mm512_add_epi32(d, g);
    let exceeds_q = _mm512_cmpge_epu32_mask(v, q);
    _mm512_mask_sub_epi32(v, exceeds_q, v, q)
}

#[inline(always)]
pub unsafe fn barrett_mul_4194304_759304193(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(759304193);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(23724876));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(4194304));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-759304193));
    let v = _mm512_add_epi32(d, g);
    let exceeds_q = _mm512_cmpge_epu32_mask(v, q);
    _mm512_mask_sub_epi32(v, exceeds_q, v, q)
}

// // (a + b*c)% (2^32-99)
// #[inline(always)]
// pub unsafe fn mla_mod32(a: __m512i, b: __m512i, c: i32) -> __m512i {
//     let mut ain = [0i32; 16];
//     let mut bin = [0i32; 16];
//     let mut res = [0i32; 16];
//     let p = (1i64 << 32) - 99;
//     _mm512_storeu_si512(ain.as_mut_ptr() as *mut __m512i, a);
//     _mm512_storeu_si512(bin.as_mut_ptr() as *mut __m512i, b);
//     for i in 0..16 {
//         let prod = (ain[i] as i64) + (bin[i] as i64) * (c as i64);
//         res[i] = (prod % p) as i32;
//     }
//     _mm512_loadu_si512(res.as_ptr() as *const __m512i)
// }

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn mla_mod32(a: __m512i, v1: __m512i) -> __m512i {
    let mut ain = [0u32; 16];
    let mut v1in = [0u32; 16];
    let mut res = [0i32; 16];
    
    let p0 = 759207937i64;
    let p1 = 759304193i64;
    let m_rns = p0 * p1;
    let half_m = m_rns / 2;
    let q = (1i64 << 32) - 99;
    
    _mm512_storeu_si512(ain.as_mut_ptr() as *mut __m512i, a);
    _mm512_storeu_si512(v1in.as_mut_ptr() as *mut __m512i, v1);
    
    for i in 0..16 {
        let fv = (ain[i] as i64) + (v1in[i] as i64) * p0;
        if fv > half_m {
            let fv_q = fv % q;
            let m_q = m_rns % q;
            res[i] = ((fv_q + q - m_q) % q) as i32;
        } else {
            res[i] = (fv % q) as i32;
        }
    }
    _mm512_loadu_si512(res.as_ptr() as *const __m512i)
}

// #[inline(always)]
// pub unsafe fn mla_mod32(a: __m512i, v1: __m512i) -> __m512i {
//     let q_val = 4294967197u32; // 2^32 - 99
//     let q = _mm512_set1_epi32(q_val as i32);
//     let p0 = _mm512_set1_epi32(759207937);

//     // ==========================================
//     // 步驟 1: 計算 fv_mod_q = (a + v1 * p0) % q
//     // ==========================================
//     let lo = _mm512_mullo_epi32(v1, p0);
//     // 這裡直接呼叫你寫好的神級 mulhi 函式
//     let hi = _mm512_mulhi_epu32(v1, p0); 
    
//     // 魔法：Modulo 2^32 - 99 等同於 hi * 99 + lo
//     let hi99 = _mm512_mullo_epi32(hi, _mm512_set1_epi32(99));
//     let term2 = _mm512_add_epi32(hi99, a); 
//     let sum = _mm512_add_epi32(lo, term2);

//     // 處理 32-bit 加法溢位 (如果 sum < term2 代表溢位，補上 99)
//     let carry = _mm512_cmpgt_epu32_mask(term2, sum);
//     let sum2 = _mm512_mask_add_epi32(sum, carry, sum, _mm512_set1_epi32(99));

//     // 如果結果 >= q，減去 q
//     let ge_q = _mm512_cmpge_epu32_mask(sum2, q);
//     let fv_mod_q = _mm512_mask_sub_epi32(sum2, ge_q, sum2, q);

//     // ==========================================
//     // 步驟 2: 判斷是否大於 half_m
//     // 條件等價於: v1 > 379652096 | (v1 == 379652096 & a >= 379603969)
//     // ==========================================
//     let gt_v1 = _mm512_cmpgt_epu32_mask(v1, _mm512_set1_epi32(379652096));
//     let eq_v1 = _mm512_cmpeq_epi32_mask(v1, _mm512_set1_epi32(379652096));
//     let ge_a = _mm512_cmpge_epu32_mask(a, _mm512_set1_epi32(379603969));
//     // is_exceed 就是 fv > half_m 的 Mask!
//     let is_exceed = gt_v1 | (eq_v1 & ge_a);

//     // ==========================================
//     // 步驟 3: 如果 is_exceed，修正數值
//     // ==========================================
//     // q - (m_rns % q) = 4294967197 - 459965471 = 3835001726
//     let to_add = _mm512_set1_epi32(3835001726u32 as i32);
//     let fv_exceed = _mm512_add_epi32(fv_mod_q, to_add);

//     // 同樣處理溢位
//     let carry_exceed = _mm512_cmpgt_epu32_mask(to_add, fv_exceed);
//     let fv_exceed2 = _mm512_mask_add_epi32(fv_exceed, carry_exceed, fv_exceed, _mm512_set1_epi32(99));

//     let ge_q2 = _mm512_cmpge_epu32_mask(fv_exceed2, q);
//     let fv_final_exceed = _mm512_mask_sub_epi32(fv_exceed2, ge_q2, fv_exceed2, q);

//     // ==========================================
//     // 步驟 4: 根據 Mask 混合最終結果
//     // ==========================================
//     _mm512_mask_blend_epi32(is_exceed, fv_mod_q, fv_final_exceed)
// }

#[cfg(test)]
mod tests {
    use super::*;
    use std::arch::x86_64::*;

    #[test]
    fn test_barrett_mul_759304193() {
        if !is_x86_feature_detected!("avx512f") {
            println!("Skipping AVX-512 test: CPU feature not detected.");
            return;
        }

        unsafe {
            let modulus: u64 = 759304193;
            let multiplier: u64 = 161561972;
            
            // 測試案例，包含邊界值與隨機大整數
            let input_data: [u32; 16] = [
                0, 1, 100, 759304192, 
                759304193, 759304194, 1000000000, 2000000000,
                3000000000, 4000000000, 4294967295, 12345678,
                87654321, 55555555, 99999999, 10101010
            ];

            // 修正點：將指標先轉為 *const i32，最後必須轉為 *const __m512i
            let a_v = _mm512_loadu_si512(input_data.as_ptr() as *const __m512i);

            // 執行 Hachi 論文中提到的優化邏輯
            let result_v = barrett_mul_759304193(a_v);

            let mut output_data = [0u32; 16];
            // 修正點：同理，存儲時必須轉為 *mut __m512i
            _mm512_storeu_si512(output_data.as_mut_ptr() as *mut __m512i, result_v);

            for i in 0..16 {
                let expected = ((input_data[i] as u64 * multiplier) % modulus) as u32;
                assert_eq!(
                    ((output_data[i] as u64) % modulus) as u32, 
                    expected, 
                    "驗證失敗！輸入: {}, 預期: {}, 實際: {}", 
                    input_data[i], expected, output_data[i]
                );
            }
        }
    }
}