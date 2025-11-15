#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;

// transpose a 16x16 matrix
pub unsafe fn transpose(a: [__m256i;16]) -> [__m256i;16]{
    let mut t0 = [_mm256_setzero_si256(); 16];
    let mut t1 = [_mm256_setzero_si256(); 16];
    let mut t2 = [_mm256_setzero_si256(); 16];
    let mut out = [_mm256_setzero_si256(); 16];

    for i in 0..8 {
        let a_even = a[i * 2];
        let a_odd = a[i * 2 + 1];
        t0[i * 2]     = _mm256_unpacklo_epi16(a_even, a_odd);
        t0[i * 2 + 1] = _mm256_unpackhi_epi16(a_even, a_odd);
    }

    for i in 0..4 {
        let base = i * 4;
        let t0_0 = t0[base];
        let t0_1 = t0[base + 1];
        let t0_2 = t0[base + 2];
        let t0_3 = t0[base + 3];
        t1[base]     = _mm256_unpacklo_epi32(t0_0, t0_2);
        t1[base + 1] = _mm256_unpackhi_epi32(t0_0, t0_2);
        t1[base + 2] = _mm256_unpacklo_epi32(t0_1, t0_3);
        t1[base + 3] = _mm256_unpackhi_epi32(t0_1, t0_3);
    }

    for i in 0..2 {
        let base = i * 8;
        let t1_0 = t1[base];
        let t1_1 = t1[base + 1];
        let t1_2 = t1[base + 2];
        let t1_3 = t1[base + 3];
        let t1_4 = t1[base + 4];
        let t1_5 = t1[base + 5];
        let t1_6 = t1[base + 6];
        let t1_7 = t1[base + 7];
        t2[base]     = _mm256_unpacklo_epi64(t1_0, t1_4);
        t2[base + 1] = _mm256_unpackhi_epi64(t1_0, t1_4);
        t2[base + 2] = _mm256_unpacklo_epi64(t1_1, t1_5);
        t2[base + 3] = _mm256_unpackhi_epi64(t1_1, t1_5);
        t2[base + 4] = _mm256_unpacklo_epi64(t1_2, t1_6);
        t2[base + 5] = _mm256_unpackhi_epi64(t1_2, t1_6);
        t2[base + 6] = _mm256_unpacklo_epi64(t1_3, t1_7);
        t2[base + 7] = _mm256_unpackhi_epi64(t1_3, t1_7);
    }

    for i in 0..8 {
        let a_lo = t2[i];
        let a_hi = t2[i + 8];
        out[i] = _mm256_permute2x128_si256(a_lo, a_hi, 0x20);
        out[i+8] = _mm256_permute2x128_si256(a_lo, a_hi, 0x31);
    }
    out
}