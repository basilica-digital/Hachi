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

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
// montgomery multiplication
pub unsafe fn montproduct_759207937(x: __m512i, y: __m512i) -> __m512i{
    let lo = _mm512_mullo_epi32(x, y);
    let hi = _mm512_mulhi_epu32(x, y);
    let d = _mm512_mullo_epi32(lo, _mm512_set1_epi32(754935809));
    let e = _mm512_mulhi_epu32(d, _mm512_set1_epi32(759207937));
    let res = _mm512_sub_epi32(hi, e);
    let is_negative = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), res);
    _mm512_mask_add_epi32(res, is_negative, res, _mm512_set1_epi32(759207937))
}

#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
#[inline(always)]
pub unsafe fn montproduct_759304193(x: __m512i, y: __m512i) -> __m512i{
    let lo = _mm512_mullo_epi32(x, y);
    let hi = _mm512_mulhi_epu32(x, y);
    let d = _mm512_mullo_epi32(lo, _mm512_set1_epi32(331214849));
    let e = _mm512_mulhi_epu32(d, _mm512_set1_epi32(759304193));
   let res = _mm512_sub_epi32(hi, e);
    let is_negative = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), res);
    _mm512_mask_add_epi32(res, is_negative, res, _mm512_set1_epi32(759304193))
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

#[inline(always)]
pub unsafe fn mla_mod32(a: __m512i, v1: __m512i) -> __m512i {
    let mut ain = [0u32; 16]; // 建議用無號讀入分量
    let mut v1in = [0u32; 16];
    let mut res = [0i32; 16];
    
    let p0 = 759207937i64;
    let p1 = 759304193i64;
    let m_rns = p0 * p1;        // RNS 總容量 (約 2^59)
    let half_m = m_rns / 2;     // 判定正負的邊界
    let q = (1i64 << 32) - 99;  // 目標模數
    
    _mm512_storeu_si512(ain.as_mut_ptr() as *mut __m512i, a);
    _mm512_storeu_si512(v1in.as_mut_ptr() as *mut __m512i, v1);
    
    for i in 0..16 {
        // 1. 還原出 [0, p0*p1) 之間的完整整數 fv
        let fv = (ain[i] as i64) + (v1in[i] as i64) * p0;
        
        // 2. 符號校正與縮減
        if fv > half_m {
            // 處理負數情況: fv 在數學上代表 (fv - m_rns)
            // 我們要計算 (fv - m_rns) % q
            // 為了避免處理負數取模的混亂，可以使用:
            // (fv % q + (q - (m_rns % q))) % q
            let fv_q = fv % q;
            let m_q = m_rns % q;
            res[i] = ((fv_q + q - m_q) % q) as i32;
        } else {
            // 處理正數情況
            res[i] = (fv % q) as i32;
        }
    }
    _mm512_loadu_si512(res.as_ptr() as *const __m512i)
}


// #[cfg(test)]
// mod tests {
//     use super::*;
//     use core::arch::x86_64::*;

//     unsafe fn dump_m512(v: __m512i) -> [u32; 16] {
//         let mut arr = [0u32; 16];
//         _mm512_storeu_si512(arr.as_mut_ptr() as *mut u32 as *mut __m512i, v);
//         arr
//     }

//     macro_rules! test_barrett {
//         ($func_name:ident, $q:expr) => {
//             #[test]
//             fn $func_name() {
//                 if !is_x86_feature_detected!("avx512f") { return; }
//                 use rand::Rng;

//                 unsafe {
//                     let q = $q as u64;
//                     let mut rng = rand::thread_rng();
//                     let mut inputs = [0u32; 16];
                    
//                     for _ in 0..1000000 {
//                         for i in 0..16 {
//                             inputs[i] = rng.gen::<u32>();
//                         }

//                         let input_vec = _mm512_loadu_si512(inputs.as_ptr() as *const __m512i);
//                         let output_vec = super::$func_name(input_vec);
//                         let outputs = dump_m512(output_vec);

//                         for i in 0..16 {
//                             assert_eq!(
//                                 inputs[i] as u64 % q, 
//                                 outputs[i] as u64 % q, 
//                                 "Fail at index {}, input: {}, q: {}", i, inputs[i], q
//                             );
//                         }
//                     }
//                 }
//             }
//         };
//     }

//     test_barrett!(barrett_fake_759207937, 759207937);
//     test_barrett!(barrett_fake_759304193, 759304193);
// }


// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::arch::x86_64::*;

//     unsafe fn set_all_epi32(val: u32) -> __m512i {
//         _mm512_set1_epi32(val as i32)
//     }

//     unsafe fn get_first_u32(vec: __m512i) -> u32 {
//         let low_128 = _mm512_extracti32x4_epi32::<0>(vec);
//         _mm_cvtsi128_si32(low_128) as u32
//     }

//     // 簡單的隨機數產生器 (Xorshift)
//     fn xorshift32(state: &mut u32) -> u32 {
//         let mut x = *state;
//         x ^= x << 13;
//         x ^= x >> 17;
//         x ^= x << 5;
//         *state = x;
//         x
//     }

//     #[test]
//     fn test_montgomery_massive_random() {
//         if is_x86_feature_detected!("avx512f") {
//             unsafe {
//                 let mut seed: u32 = 0x12345678;
//                 let num_tests = 1_000_000;
                
//                 let q1: u32 = 759207937;
//                 let q2: u32 = 759304193;
//                 let r_mod: u128 = (1u128 << 32);

//                 println!("開始進行 {} 組隨機測資驗證...", num_tests);

//                 for _ in 0..num_tests {
//                     let x = xorshift32(&mut seed) % q1;
//                     let y = xorshift32(&mut seed) % q1;

//                     // --- 測試 759207937 ---
//                     let res_vec1 = montproduct_759207937(set_all_epi32(x), set_all_epi32(y));
//                     let raw_res1 = get_first_u32(res_vec1);
                    
//                     let res1 = if (raw_res1 as i32) < 0 {
//                         raw_res1.wrapping_add(q1)
//                     } else {
//                         raw_res1
//                     } as u128;

//                     let left1 = (res1 * r_mod) % q1 as u128;
//                     let right1 = (x as u128 * y as u128) % q1 as u128;
                    
//                     if left1 != right1 {
//                         panic!("759207937 失敗! x: {}, y: {}, res: {}, left: {}, right: {}", x, y, res1, left1, right1);
//                     }

//                     // --- 測試 759304193 ---
//                     let res_vec2 = montproduct_759304193(set_all_epi32(x), set_all_epi32(y));
//                     let raw_res2 = get_first_u32(res_vec2);
                    
//                     let res2 = if (raw_res2 as i32) < 0 {
//                         raw_res2.wrapping_add(q2)
//                     } else {
//                         raw_res2
//                     } as u128;

//                     let left2 = (res2 * r_mod) % q2 as u128;
//                     let right2 = (x as u128 * y as u128) % q2 as u128;

//                     if left2 != right2 {
//                         panic!("759304193 失敗! x: {}, y: {}, res: {}, left: {}, right: {}", x, y, res2, left2, right2);
//                     }
//                 }

//                 println!("✅ 通過所有 {} 組隨機測資驗證！", num_tests);
//             }
//         } else {
//             println!("跳過測試：CPU 不支援 AVX-512");
//         }
//     }
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