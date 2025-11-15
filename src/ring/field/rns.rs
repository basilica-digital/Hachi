#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
use crate::NTT::f26113::*;
use crate::NTT::f25601::*;
use crate::NTT::f23041::*;
use crate::NTT::f19457::*;
use crate::NTT::f18433::*;
//[26113, 25601, 23041, 19457, 18433]

// 32x8 -> 5x 16x16
pub unsafe fn rns_ntt(a: [__m256i;32]) -> [__m256i;32]{
    let mut poly_16: [[__m256i;16];5] = [[_mm256_setzero_si256(); 16];5];
    for i in 0..16{
        poly_16[0][i] = reduce_26113(a[2*i], a[2*i+1]);
        poly_16[1][i] = reduce_25601(a[2*i], a[2*i+1]);
        poly_16[2][i] = reduce_23041(a[2*i], a[2*i+1]);
        poly_16[3][i] = reduce_19457(a[2*i], a[2*i+1]);
        poly_16[4][i] = reduce_18433(a[2*i], a[2*i+1]);
    }
    poly_16[0] = ntt_26113(poly_16[0]);
    poly_16[1] = ntt_25601(poly_16[1]);
    poly_16[2] = ntt_23041(poly_16[2]);
    poly_16[3] = ntt_19457(poly_16[3]);
    poly_16[4] = ntt_18433(poly_16[4]);
}

// reduce
pub unsafe fn reduce_26113(a: __m256i, b: __m256i) -> __m256i {
    let a1 = barrett_fake_26113_32(a);
    let b1 = barrett_fake_26113_32(b);
    let c = _mm256_packs_epi32(a, b);
    let res = _mm256_permute4x64_epi64(c, 0xD8);
    res
}

pub unsafe fn reduce_25601(a: __m256i, b: __m256i) -> __m256i {
    let a1 = barrett_fake_25601_32(a);
    let b1 = barrett_fake_25601_32(b);
    let c = _mm256_packs_epi32(a, b);
    let res = _mm256_permute4x64_epi64(c, 0xD8);
    res
}

pub unsafe fn reduce_23041(a: __m256i, b: __m256i) -> __m256i {
    let a1 = barrett_fake_23041_32(a);
    let b1 = barrett_fake_23041_32(b);
    let c = _mm256_packs_epi32(a, b);
    let res = _mm256_permute4x64_epi64(c, 0xD8);
    res
}

pub unsafe fn reduce_19457(a: __m256i, b: __m256i) -> __m256i {
    let a1 = barrett_fake_19457_32(a);
    let b1 = barrett_fake_19457_32(b);
    let c = _mm256_packs_epi32(a, b);
    let res = _mm256_permute4x64_epi64(c, 0xD8);
    res
}

pub unsafe fn reduce_18433(a: __m256i, b: __m256i) -> __m256i {
    let a1 = barrett_fake_18433_32(a);
    let b1 = barrett_fake_18433_32(b);
    let c = _mm256_packs_epi32(a, b);
    let res = _mm256_permute4x64_epi64(c, 0xD8);
    res
}