#[inline(never)]
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
    let diff = _mm512_sub_epi32(b, a);
	let v = barrett_mul_759304193(diff, _mm512_set1_epi32(161561972));
    mla_mod32(a, v, 759207937)
}

pub unsafe fn irns(val1: &[__m512i; 64], val2: &[__m512i; 64], out: *mut __m512i, offset: usize){
    for i in 0..64{
        _mm512_storeu_si512(out.add(offset+i) as *mut _, irns_vec(val1[i], val2[i]));
    }
}