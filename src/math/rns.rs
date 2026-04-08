use std::arch::x86_64::*;
use crate::math::field_simd::*;

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

#[inline(always)]
pub unsafe fn rns_2048(input_data: &[u32], d1_target: &mut [u32], d2_target: &mut [u32]) {
    let x = load_vec_to_m512(input_data, 0);
    let res1 = barrett_fake_2079301633(x);
    let res2 = barrett_fake_2079305729(x);
    store_m512_to_vec(res1, d1_target, 0);
    store_m512_to_vec(res2, d2_target, 0);
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

pub unsafe fn irns_vec_2048(a: __m512i, b: __m512i) -> __m512i {
    let mut diff = _mm512_sub_epi32(b, a);
    let is_neg = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), diff);
    diff = _mm512_mask_add_epi32(diff, is_neg, diff, _mm512_set1_epi32(2079305729));
    let v = barrett_mul_2079305729(diff);
    mla_mod32_2048(a, v)
}

pub unsafe fn irns_2048(val1: &[__m512i; 128], val2: &[__m512i; 128], out: *mut __m512i, offset: usize) {
    for i in 0..128 {
        _mm512_storeu_si512(out.add(offset + i) as *mut _, irns_vec_2048(val1[i], val2[i]));
    }
}