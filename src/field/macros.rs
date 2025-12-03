#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;

#[target_feature(enable = "avx2")]
pub unsafe fn _mm256_mulhi_epi32(a: __m256i, b: __m256i) -> __m256i{
    let prod02 = _mm256_mul_epu32(a, b);
    let prod13 = _mm256_mul_epu32(_mm256_shuffle_epi32(a, 0xf5),_mm256_shuffle_epi32(b, 0xf5));
    let ans = _mm256_unpackhi_epi64(_mm256_unpacklo_epi32(prod02, prod13),_mm256_unpackhi_epi32(prod02, prod13));
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn pack_32_16(a: __m256i, b: __m256i) -> __m256i {
    let mask = _mm256_setr_epi8(
        0, 1, 4, 5, 8, 9, 12, 13,
        -1, -1, -1, -1, -1, -1, -1, -1,
        0, 1, 4, 5, 8, 9, 12, 13,
        -1, -1, -1, -1, -1, -1, -1, -1
    );
    let sa = _mm256_shuffle_epi8(a, mask);
    let sb = _mm256_shuffle_epi8(b, mask);
    let mixed = _mm256_unpacklo_epi64(sa, sb);
    let ret = _mm256_permute4x64_epi64(mixed, 0xD8);
    
    ret
}