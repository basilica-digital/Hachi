#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
// RNS: [7681, 10753, 11777, 12289, 13313, 15361]

#[target_feature(enable = "avx2")]
pub unsafe fn rns_ntt(a: [__m256i;32]) -> [__m256i;32]{
    let mut poly_16: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    for i in 0..16{
        poly_16[0][i] = barrett_fake_257_32(a[2*i], a[2*i+1]);
        poly_16[1][i] = barrett_fake_3329_32(a[2*i], a[2*i+1]);
        poly_16[2][i] = barrett_fake_7681_32(a[2*i], a[2*i+1]);
        poly_16[3][i] = barrett_fake_7937_32(a[2*i], a[2*i+1]);
        poly_16[4][i] = barrett_fake_9473_32(a[2*i], a[2*i+1]);
        poly_16[5][i] = barrett_fake_10753_32(a[2*i], a[2*i+1]);
    }
    poly_16[0] = ntt_26113(poly_16[0]);
    poly_16[1] = ntt_25601(poly_16[1]);
    poly_16[2] = ntt_23041(poly_16[2]);
    poly_16[3] = ntt_19457(poly_16[3]);
    poly_16[4] = ntt_18433(poly_16[4]);
    poly_16[5] = ntt_18433(poly_16[5]);
}

#[target_feature(enable = "avx2")]
pub unsafe fn rns_dotproduct(a: [[__m256i;16];6]) -> [[__m256i;16];6]{
    let mut poly_16: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    poly_16[0] = dot_product_26113(a[0]);
    poly_16[1] = dot_product_25601(a[1]);
    poly_16[2] = dot_product_23041(a[2]);
    poly_16[3] = dot_product_19457(a[3]);
    poly_16[4] = dot_product_18433(a[4]);
    poly_16[5] = dot_product_18433(a[5]);
    poly_16
}

#[target_feature(enable = "avx2")]
pub unsafe fn rns_intt(a: [__m256i;32]) -> [__m256i;32]{
    let mut v: [[__m256i;16];6] = [[_mm256_setzero_si256(); 16];6];
    let f: [i16;6] = [7681, 10753, 11777, 12289, 13313, 15361];
    for i in 0..6{
        for j in 0..16{
            v[i][j] = a[i][j];
        }
    }
    for k in 0..16{
        v[1][k] = barrett_mul_10753(_mm256_sub_epi16(v[1][k], v[0]), _mm256_set1_epi16(5380));
        v[2][k] = barrett_mul_11777(_mm256_sub_epi16(v[2][k], v[0]), _mm256_set1_epi16(1475));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[0]), _mm256_set1_epi16(4099));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[0]), _mm256_set1_epi16(7264));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[0]), _mm256_set1_epi16(2));
        v[2][k] = barrett_mul_11777(_mm256_sub_epi16(v[2][k], v[1]), _mm256_set1_epi16(5900));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[1]), _mm256_set1_epi16(8));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[1]), _mm256_set1_epi16(7993));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[1]), _mm256_set1_epi16(10244));
        v[3][k] = barrett_mul_12289(_mm256_sub_epi16(v[3][k], v[2]), _mm256_set1_epi16(24));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[2]), _mm256_set1_epi16(8884));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[2]), _mm256_set1_epi16(8782));
        v[4][k] = barrett_mul_13313(_mm256_sub_epi16(v[4][k], v[3]), _mm256_set1_epi16(13));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[3]), _mm256_set1_epi16(5));
        v[5][k] = barrett_mul_15361(_mm256_sub_epi16(v[5][k], v[4]), _mm256_set1_epi16(7688));
    }
    let mut u: [[__m256i; 32]; 6] = [[_mm256_setzero_si256(); 32];6];
    for i in 0..6{
        for j in 0..16{
            v[i][j] = _mm256_add_epi16(v[i][j], _mm256_set1_epi16(f[i]));
            u[i][2*j] = _mm256_cvtepu16_epi32(_mm256_castsi256_si128(v[i][j]));
            u[i][2*j+1] = _mm256_cvtepu16_epi32(_mm256_extracti128_si256::<1>(v[i][j]));
        }
    }
    let w: [i32;6] = [1, 7681, 82593793, 2044513639, -529992779, 837237844, 1678733866];
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