#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::fields::*;
// RNS: [257, 3329, 7681, 7937, 9473, 10753]

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