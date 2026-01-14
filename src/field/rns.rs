#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::ntt::{ntt::*, intt::*, dotproduct::*};
// RNS: [7681, 10753, 11777, 12289, 13313, 15361]

#[target_feature(enable = "avx2")]
pub unsafe fn rns_ntt(a: [__m256i;32]) -> [[__m256i; 16]; 6]{
    let mut poly_16: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    for i in 0..16{
        poly_16[0][i] = barrett_fake_7681_32(a[2*i], a[2*i+1]);
        poly_16[1][i] = barrett_fake_10753_32(a[2*i], a[2*i+1]);
        poly_16[2][i] = barrett_fake_11777_32(a[2*i], a[2*i+1]);
        poly_16[3][i] = barrett_fake_12289_32(a[2*i], a[2*i+1]);
        poly_16[4][i] = barrett_fake_13313_32(a[2*i], a[2*i+1]);
        poly_16[5][i] = barrett_fake_15361_32(a[2*i], a[2*i+1]);
    }
    poly_16[0] = ntt_7681(poly_16[0]);
    poly_16[1] = ntt_10753(poly_16[1]);
    poly_16[2] = ntt_11777(poly_16[2]);
    poly_16[3] = ntt_12289(poly_16[3]);
    poly_16[4] = ntt_13313(poly_16[4]);
    poly_16[5] = ntt_15361(poly_16[5]);
    poly_16
}

#[target_feature(enable = "avx2")]
pub unsafe fn rns_dotproduct(a: [[__m256i;16];6], b: [[__m256i;16];6]) -> [[__m256i;16];6]{
    let mut poly_16: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    poly_16[0] = dot_product_7681(a[0], b[0]);
    poly_16[1] = dot_product_10753(a[1], b[1]);
    poly_16[2] = dot_product_11777(a[2], b[2]);
    poly_16[3] = dot_product_12289(a[3], b[3]);
    poly_16[4] = dot_product_13313(a[4], b[4]);
    poly_16[5] = dot_product_15361(a[5], b[5]);
    poly_16
}

#[target_feature(enable = "avx2")]
pub unsafe fn rns_intt(a: [[__m256i; 16]; 6]) -> [__m256i;32]{
    let mut v: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    let f: [i16;6] = [7681, 10753, 11777, 12289, 13313, 15361];
    for i in 0..6{
        for j in 0..16{
            v[i][j] = a[i][j];
        }
    }
    for k in 0..16{
        v[1][k] = barrett_mul_10753(_mm256_sub_epi16(v[1][k], v[0][k]), _mm256_set1_epi16(-5373), _mm256_set1_epi16(-32747));
        v[2][k] = barrett_mul_11777(_mm256_sub_epi16(v[2][k], v[0][k]), _mm256_set1_epi16(1475), _mm256_set1_epi16(8208));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[0][k]), _mm256_set1_epi16(4099), _mm256_set1_epi16(21860));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[0][k]), _mm256_set1_epi16(-6049), _mm256_set1_epi16(-29777));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[0][k]), _mm256_set1_epi16(2), _mm256_set1_epi16(9));
        v[2][k] = barrett_mul_11777(_mm256_sub_epi16(v[2][k], v[1][k]), _mm256_set1_epi16(-5877), _mm256_set1_epi16(-32704));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[1][k]), _mm256_set1_epi16(8), _mm256_set1_epi16(43));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[1][k]), _mm256_set1_epi16(-5320), _mm256_set1_epi16(-26189));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[1][k]), _mm256_set1_epi16(-5117), _mm256_set1_epi16(-21831));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[2][k]), _mm256_set1_epi16(24), _mm256_set1_epi16(128));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[2][k]), _mm256_set1_epi16(-4429), _mm256_set1_epi16(-21803));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[2][k]), _mm256_set1_epi16(-6579), _mm256_set1_epi16(-28069));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[3][k]), _mm256_set1_epi16(13), _mm256_set1_epi16(64));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[3][k]), _mm256_set1_epi16(5), _mm256_set1_epi16(21));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[4][k]), _mm256_set1_epi16(-7673), _mm256_set1_epi16(-32736));
    }
    let mut u: [[__m256i; 32]; 6] = [[_mm256_setzero_si256(); 32];6];
    for i in 0..6{
        for j in 0..16{
            v[i][j] = _mm256_add_epi16(v[i][j], _mm256_set1_epi16(f[i]));
            u[i][2*j] = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(v[i][j]));
            u[i][2*j+1] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256::<1>(v[i][j]));
        }
    }
    let w: [i32;7] = [1, 7681, 82593793, 2044513639, -529992779, 837237844, 1678733866];
    let mut res: [[__m256i;32];2] = [[_mm256_setzero_si256();32];2];
    for i in 0..6{
        for j in 0..32{
            let tmp = _mm256_mul_epi32(v[i][j], w[i][j]);
            res[i][0] = _mm256_add_epi64(res[i][0], tmp);
            let v_shuf = _mm256_shuffle_epi32::<0xB1>(v[i][j]);
            let w_shuf = _mm256_shuffle_epi32::<0xB1>(w[i][j]);
            let tmp = _mm256_mul_epi32(v_shuf, w_shuf);
            res[i][1] = _mm256_add_epi64(res[i][1], tmp);
        }
    }
    // reduction
    let mut ans: [__m256i;32] = [_mm256_setzero_si256();32];
    for i in 0..32{
        ans[i] = reduce32(res[i][0], res[i][1]);
    }
    ans
}


// #[cfg(test)]
// mod tests {
//     use super::*;
//     use std::arch::x86_64::*;

//     // 輔助工具：將 __m256i 陣列轉回 i32 陣列以便比對
//     unsafe fn dump_to_i32_256(data: [__m256i; 32]) -> [i32; 256] {
//         let mut res = [0i32; 256];
//         for i in 0..32 {
//             let ptr = &data[i] as *const __m256i as *const i32;
//             for j in 0..8 {
//                 res[i * 8 + j] = *ptr.add(j);
//             }
//         }
//         res
//     }

//     // 輔助工具：將 i32 數據填入 __m256i 陣列
//     unsafe fn load_from_i32_256(data: &[i32; 256]) -> [__m256i; 32] {
//         let mut res = [_mm256_setzero_si256(); 32];
//         for i in 0..32 {
//             res[i] = _mm256_loadu_si256(data[i * 8..].as_ptr() as *const __m256i);
//         }
//         res
//     }

//     #[test]
//     fn test_rns_ntt_intt_identity() {
//         if !is_x86_feature_detected!("avx2") { return; }

//         unsafe {
//             let step = 256;
//             let mut base_val: i64 = i16::MIN as i64; 

//             while base_val <= i16::MAX as i64 {
//                 let mut input_raw = [0i32; 256];
//                 for i in 0..256 {
//                     input_raw[i] = ((base_val + i as i64) % 3329) as i32; 
//                 }

//                 let a = load_from_i32_256(&input_raw);

//                 let forward_rns = rns_ntt(a);
//                 let backward_a = rns_intt(forward_rns);

//                 let output_raw = dump_to_i32_256(backward_a);

//                 for i in 0..256 {
//                     assert_eq!(
//                         input_raw[i], 
//                         output_raw[i], 
//                         "\nIdentity Failed!\nIndex: {}\nExpected: {}\nGot: {}\nBase: {}", 
//                         i, input_raw[i], output_raw[i], base_val
//                     );
//                 }

//                 base_val += step as i64;
//             }
//             println!("SUCCESS: RNS NTT/INTT identity verified for all sampled ranges.");
//         }
//     }
// }
