use std::arch::x86_64::*;

pub type EF8 = (__m512i, __m512i, __m512i, __m512i);

#[inline(always)]
pub unsafe fn load_vec_to_m512(src: &[u32], offset: usize) -> __m512i {
    let ptr = src.as_ptr().add(offset) as *const __m512i;
    _mm512_load_si512(ptr)
}

#[inline(always)]
pub unsafe fn store_m512_to_vec(src: __m512i, dst: &mut [u32], offset: usize) {
    let ptr = dst.as_mut_ptr().add(offset) as *mut __m512i;
    _mm512_store_si512(ptr, src)
}

#[inline(always)]
pub unsafe fn _mm512_mulhi_epu32(a: __m512i, b: __m512i) -> __m512i {
    let mul_evens = _mm512_mul_epu32(a, b);
    let a_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(a);
    let b_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(b);
    let mul_odds = _mm512_mul_epu32(a_odds, b_odds);
    let hi_evens = _mm512_srli_epi64::<32>(mul_evens);
    let mask = 0xAAAA; 
    _mm512_mask_blend_epi32(mask, hi_evens, mul_odds)
}

#[inline(always)]
pub unsafe fn _mm512_mulhi_epu32_const(a: __m512i, b: __m512i) -> __m512i {
    let mul_evens = _mm512_mul_epu32(a, b);
    let a_odds = _mm512_shuffle_epi32::<0b11_11_01_01>(a);
    let mul_odds = _mm512_mul_epu32(a_odds, b);
    let hi_evens = _mm512_srli_epi64::<32>(mul_evens);
    let mask = 0xAAAA; 
    _mm512_mask_blend_epi32(mask, hi_evens, mul_odds)
}

// barrett_fake
#[inline(always)]
pub unsafe fn barrett_fake_759207937(x: __m512i) -> __m512i{
    let vec = _mm512_set1_epi32(759207937);
    let d = _mm512_mulhi_epu32(x, _mm512_set1_epi32(6));
    let e = _mm512_mullo_epi32(d, vec);
    _mm512_sub_epi32(_mm512_add_epi32(x, vec), e)
}


#[inline(always)]
pub unsafe fn barrett_fake_759304193(x: __m512i) -> __m512i{
    let vec = _mm512_set1_epi32(759304193);
    let d = _mm512_mulhi_epu32(x, _mm512_set1_epi32(6));
    let e = _mm512_mullo_epi32(d, vec);
    _mm512_sub_epi32(_mm512_add_epi32(x, vec), e)
}


#[inline(always)]
pub unsafe fn barrett_fake_2079301633(x: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(2079301633);
    let ge_q1 = _mm512_cmpge_epu32_mask(x, q);
    let x1 = _mm512_mask_sub_epi32(x, ge_q1, x, q);
    let ge_q2 = _mm512_cmpge_epu32_mask(x1, q);
    _mm512_mask_sub_epi32(x1, ge_q2, x1, q)
}


#[inline(always)]
pub unsafe fn barrett_fake_2079305729(x: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(2079305729);
    let ge_q1 = _mm512_cmpge_epu32_mask(x, q);
    let x1 = _mm512_mask_sub_epi32(x, ge_q1, x, q);
    let ge_q2 = _mm512_cmpge_epu32_mask(x1, q);
    _mm512_mask_sub_epi32(x1, ge_q2, x1, q)
}


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


#[inline(always)]
pub unsafe fn montproduct(x: __m512i, y: __m512i, q: __m512i, q_inv: __m512i) -> __m512i{
    let lo = _mm512_mullo_epi32(x, y);
    let hi = _mm512_mulhi_epu32(x, y);
    let d = _mm512_mullo_epi32(lo, q_inv);
    let e = _mm512_mulhi_epu32_const(d, q);
   let res = _mm512_sub_epi32(hi, e);
    let is_negative = _mm512_cmpgt_epi32_mask(_mm512_setzero_si512(), res);
    _mm512_mask_add_epi32(res, is_negative, res, q)
}


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
pub unsafe fn barrett_mul_2079305729(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(2079305729);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(1048575));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(507643));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-2079305729));
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

#[inline(always)]
pub unsafe fn barrett_mul_2097152_2079301633(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(2079301633);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(4331838));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(2097152));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-2079301633));
    let v = _mm512_add_epi32(d, g);
    let exceeds_q = _mm512_cmpge_epu32_mask(v, q);
    _mm512_mask_sub_epi32(v, exceeds_q, v, q)
}

#[inline(always)]
pub unsafe fn barrett_mul_2097152_2079305729(a: __m512i) -> __m512i{
    let q = _mm512_set1_epi32(2079305729);
    let t = _mm512_mulhi_epu32(a, _mm512_set1_epi32(4331830));
    let d = _mm512_mullo_epi32(a, _mm512_set1_epi32(2097152));
    let g = _mm512_mullo_epi32(t, _mm512_set1_epi32(-2079305729));
    let v = _mm512_add_epi32(d, g);
    let exceeds_q = _mm512_cmpge_epu32_mask(v, q);
    _mm512_mask_sub_epi32(v, exceeds_q, v, q)
}


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


#[inline(always)]
pub unsafe fn mla_mod32_2048(a: __m512i, v1: __m512i) -> __m512i {
    let mut ain = [0u32; 16];
    let mut v1in = [0u32; 16];
    let mut res = [0i32; 16];
    
    let p0 = 2079301633i64;
    let p1 = 2079305729i64;
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


#[inline(always)]
pub unsafe fn mersenne_reduce(x: __m512i) -> __m512i {
    let c99 = _mm512_set1_epi64(99);
    let mask32 = _mm512_set1_epi64(0xFFFFFFFF);
    let q_vec = _mm512_set1_epi64(4294967197);

    let hi1 = _mm512_srli_epi64(x, 32);
    let lo1 = _mm512_and_si512(x, mask32);
    let hi1_99 = _mm512_mullo_epi64(hi1, c99);
    let x1 = _mm512_add_epi64(lo1, hi1_99);

    let hi2 = _mm512_srli_epi64(x1, 32);
    let lo2 = _mm512_and_si512(x1, mask32);
    let hi2_99 = _mm512_mullo_epi64(hi2, c99);
    let x2 = _mm512_add_epi64(lo2, hi2_99);

    let ge_q = _mm512_cmpge_epu64_mask(x2, q_vec);
    _mm512_mask_sub_epi64(x2, ge_q, x2, q_vec)
}

#[inline(always)]
pub unsafe fn extmul_8x(
    a0: __m512i, a1: __m512i, a2: __m512i, a3: __m512i,
    b0: __m512i, b1: __m512i, b2: __m512i, b3: __m512i
) -> (__m512i, __m512i, __m512i, __m512i) {
    let t0 = _mm512_mul_epu32(a0, b0);
    let t1 = _mm512_add_epi64(_mm512_mul_epu32(a0, b1), _mm512_mul_epu32(a1, b0));
    let t2 = _mm512_add_epi64(_mm512_add_epi64(_mm512_mul_epu32(a0, b2), _mm512_mul_epu32(a1, b1)), _mm512_mul_epu32(a2, b0));
    let t3 = _mm512_add_epi64(_mm512_add_epi64(_mm512_add_epi64(_mm512_mul_epu32(a0, b3), _mm512_mul_epu32(a1, b2)), _mm512_mul_epu32(a2, b1)), _mm512_mul_epu32(a3, b0));
    let t4 = _mm512_add_epi64(_mm512_add_epi64(_mm512_mul_epu32(a1, b3), _mm512_mul_epu32(a2, b2)), _mm512_mul_epu32(a3, b1));
    let t5 = _mm512_add_epi64(_mm512_mul_epu32(a2, b3), _mm512_mul_epu32(a3, b2));
    let t6 = _mm512_mul_epu32(a3, b3);
    let c0_raw = _mm512_add_epi64(t0, _mm512_slli_epi64(t4, 1));
    let c1_raw = _mm512_add_epi64(t1, _mm512_slli_epi64(t5, 1));
    let c2_raw = _mm512_add_epi64(t2, _mm512_slli_epi64(t6, 1));
    (
        mersenne_reduce(c0_raw),
        mersenne_reduce(c1_raw),
        mersenne_reduce(c2_raw),
        mersenne_reduce(t3)
    )
}

#[inline(always)]
pub unsafe fn accum_high_part_simd(
    dst: *mut u32, src: *const u32, idx_raw: u32, 
    lut_vec: __m512i, q_vec: __m512i, c99_vec: __m512i,
) {
    let is_pos = idx_raw < 1024;
    let k = if is_pos { idx_raw as usize } else { (idx_raw - 1024) as usize };
    if k == 0 { return; }
    
    let src_ptr = src.add(1024 - k);
    let mut i = 0;
    
    while i + 16 <= k {
        let b_vec = _mm512_loadu_si512(src_ptr.add(i) as *const _);
        let prod = _mm512_permutexvar_epi32(b_vec, lut_vec);
        let r_vec = _mm512_loadu_si512(dst.add(i) as *const _);
        
        if is_pos {
            let sum = _mm512_add_epi32(r_vec, prod);
            let carry = _mm512_cmplt_epu32_mask(sum, r_vec);
            let res = _mm512_mask_add_epi32(sum, carry, sum, c99_vec);
            let ge_q = _mm512_cmpge_epu32_mask(res, q_vec);
            let res_mod = _mm512_mask_sub_epi32(res, ge_q, res, q_vec);
            _mm512_storeu_si512(dst.add(i) as *mut _, res_mod);
        } else {
            let lt = _mm512_cmplt_epu32_mask(r_vec, prod);
            let diff = _mm512_sub_epi32(r_vec, prod);
            let res = _mm512_mask_add_epi32(diff, lt, diff, q_vec);
            _mm512_storeu_si512(dst.add(i) as *mut _, res);
        }
        i += 16;
    }
    
    if i < k {
        let tail = k - i;
        let mask = (1u16 << tail) - 1;
        
        let b_vec = _mm512_maskz_loadu_epi32(mask, src_ptr.add(i) as *const _);
        let prod = _mm512_permutexvar_epi32(b_vec, lut_vec);
        let r_vec = _mm512_maskz_loadu_epi32(mask, dst.add(i) as *const _);
        
        if is_pos {
            let sum = _mm512_add_epi32(r_vec, prod);
            let carry = _mm512_cmplt_epu32_mask(sum, r_vec);
            let res = _mm512_mask_add_epi32(sum, carry, sum, c99_vec);
            let ge_q = _mm512_cmpge_epu32_mask(res, q_vec);
            let res_mod = _mm512_mask_sub_epi32(res, ge_q, res, q_vec);
            _mm512_mask_storeu_epi32(dst.add(i) as *mut _, mask, res_mod);
        } else {
            let lt = _mm512_cmplt_epu32_mask(r_vec, prod);
            let diff = _mm512_sub_epi32(r_vec, prod);
            let res = _mm512_mask_add_epi32(diff, lt, diff, q_vec);
            _mm512_mask_storeu_epi32(dst.add(i) as *mut _, mask, res);
        }
    }
}

#[inline(always)]
pub unsafe fn u642u32(ptr: *mut [u32; 4], v0: __m512i, v1: __m512i, v2: __m512i, v3: __m512i) {
    let s0 = _mm512_cvtepi64_epi32(v0);
    let s1 = _mm512_cvtepi64_epi32(v1);
    let s2 = _mm512_cvtepi64_epi32(v2);
    let s3 = _mm512_cvtepi64_epi32(v3);

    let t0 = _mm256_unpacklo_epi32(s0, s1); 
    let t1 = _mm256_unpackhi_epi32(s0, s1); 
    let t2 = _mm256_unpacklo_epi32(s2, s3); 
    let t3 = _mm256_unpackhi_epi32(s2, s3); 

    let u0 = _mm256_unpacklo_epi64(t0, t2); 
    let u1 = _mm256_unpackhi_epi64(t0, t2); 
    let u2 = _mm256_unpacklo_epi64(t1, t3); 
    let u3 = _mm256_unpackhi_epi64(t1, t3); 

    let v01_lo = _mm256_permute2x128_si256(u0, u1, 0x20); 
    let v01_hi = _mm256_permute2x128_si256(u0, u1, 0x31); 
    let v23_lo = _mm256_permute2x128_si256(u2, u3, 0x20); 
    let v23_hi = _mm256_permute2x128_si256(u2, u3, 0x31); 

    _mm256_storeu_si256(ptr.add(0) as *mut _, v01_lo);
    _mm256_storeu_si256(ptr.add(2) as *mut _, v23_lo);
    _mm256_storeu_si256(ptr.add(4) as *mut _, v01_hi);
    _mm256_storeu_si256(ptr.add(6) as *mut _, v23_hi);
}

#[inline(always)]
pub unsafe fn base_add_8x(a: __m512i, b: __m512i) -> __m512i {
    let q = _mm512_set1_epi64(4294967197);
    let s = _mm512_add_epi64(a, b);
    let m = _mm512_cmpge_epu64_mask(s, q);
    _mm512_mask_sub_epi64(s, m, s, q)
}

#[inline(always)]
pub unsafe fn base_sub_8x(a: __m512i, b: __m512i) -> __m512i {
    let q = _mm512_set1_epi64(4294967197);
    let a_plus_q = _mm512_add_epi64(a, q);
    let diff = _mm512_sub_epi64(a_plus_q, b);
    let m = _mm512_cmpge_epu64_mask(diff, q);
    _mm512_mask_sub_epi64(diff, m, diff, q)
}

#[inline(always)]
pub unsafe fn base_mul_8x(a: __m512i, b: __m512i) -> __m512i {
    let mask32 = _mm512_set1_epi64(0xFFFFFFFF);
    let c99 = _mm512_set1_epi64(99);
    
    let prod = _mm512_mul_epu32(a, b);
    
    // Mersenne Reduce
    let hi = _mm512_srli_epi64(prod, 32);
    let lo = _mm512_and_si512(prod, mask32);
    let mut res = _mm512_add_epi64(lo, _mm512_mullo_epi64(hi, c99));
    
    // Second Reduce
    let hi2 = _mm512_srli_epi64(res, 32);
    let lo2 = _mm512_and_si512(res, mask32);
    res = _mm512_add_epi64(lo2, _mm512_mullo_epi64(hi2, c99));

    let q = _mm512_set1_epi64(4294967197);
    let ge_q = _mm512_cmpge_epu64_mask(res, q);
    _mm512_mask_sub_epi64(res, ge_q, res, q)
}

#[inline(always)]
pub unsafe fn base_sub_const_8x(a: __m512i, k: u64) -> __m512i {
    let q = _mm512_set1_epi64(4294967197);
    let k_vec = _mm512_set1_epi64(k as i64);
    let a_plus_q = _mm512_add_epi64(a, q);
    let diff = _mm512_sub_epi64(a_plus_q, k_vec);
    let m = _mm512_cmpge_epu64_mask(diff, q);
    _mm512_mask_sub_epi64(diff, m, diff, q)
}

#[inline(always)]
pub unsafe fn ext_mul_base_8x(a: EF8, b_base: __m512i) -> EF8 {
    (
        base_mul_8x(a.0, b_base),
        base_mul_8x(a.1, b_base),
        base_mul_8x(a.2, b_base),
        base_mul_8x(a.3, b_base)
    )
}

#[inline(always)]
pub unsafe fn broadcast_aos_to_soa(a: &[u32; 4]) -> EF8 {
    (
        _mm512_set1_epi64(a[0] as i64), _mm512_set1_epi64(a[1] as i64),
        _mm512_set1_epi64(a[2] as i64), _mm512_set1_epi64(a[3] as i64)
    )
}

#[inline(always)]
pub unsafe fn load8_soa(ptr: *const [u32; 4]) -> EF8 {
    let p = ptr as *const i32;
    let v0 = _mm256_setr_epi32(*p.add(0), *p.add(4), *p.add(8), *p.add(12), *p.add(16), *p.add(20), *p.add(24), *p.add(28));
    let v1 = _mm256_setr_epi32(*p.add(1), *p.add(5), *p.add(9), *p.add(13), *p.add(17), *p.add(21), *p.add(25), *p.add(29));
    let v2 = _mm256_setr_epi32(*p.add(2), *p.add(6), *p.add(10), *p.add(14), *p.add(18), *p.add(22), *p.add(26), *p.add(30));
    let v3 = _mm256_setr_epi32(*p.add(3), *p.add(7), *p.add(11), *p.add(15), *p.add(19), *p.add(23), *p.add(27), *p.add(31));
    (_mm512_cvtepu32_epi64(v0), _mm512_cvtepu32_epi64(v1), _mm512_cvtepu32_epi64(v2), _mm512_cvtepu32_epi64(v3))
}

#[inline(always)]
pub unsafe fn store8_soa(ptr: *mut [u32; 4], soa: EF8) {
    let p = ptr as *mut u32;
    let mut a0 = [0u64; 8]; let mut a1 = [0u64; 8]; let mut a2 = [0u64; 8]; let mut a3 = [0u64; 8];
    _mm512_storeu_si512(a0.as_mut_ptr() as *mut _, soa.0); _mm512_storeu_si512(a1.as_mut_ptr() as *mut _, soa.1);
    _mm512_storeu_si512(a2.as_mut_ptr() as *mut _, soa.2); _mm512_storeu_si512(a3.as_mut_ptr() as *mut _, soa.3);
    for k in 0..8 {
        *p.add(k*4 + 0) = a0[k] as u32; *p.add(k*4 + 1) = a1[k] as u32;
        *p.add(k*4 + 2) = a2[k] as u32; *p.add(k*4 + 3) = a3[k] as u32;
    }
}

#[inline(always)]
pub unsafe fn extadd_8x(a: EF8, b: EF8) -> EF8 {
    let q = _mm512_set1_epi64(4294967197);
    let add_lane = |x, y| {
        let s = _mm512_add_epi64(x, y);
        let m = _mm512_cmpge_epu64_mask(s, q);
        _mm512_mask_sub_epi64(s, m, s, q)
    };
    (add_lane(a.0, b.0), add_lane(a.1, b.1), add_lane(a.2, b.2), add_lane(a.3, b.3))
}

#[inline(always)]
pub unsafe fn extsub_8x(a: EF8, b: EF8) -> EF8 {
    let q = _mm512_set1_epi64(4294967197);
    let sub_lane = |x, y| {
        let x_plus_q = _mm512_add_epi64(x, q);
        let diff = _mm512_sub_epi64(x_plus_q, y);
        let m = _mm512_cmpge_epu64_mask(diff, q);
        _mm512_mask_sub_epi64(diff, m, diff, q)
    };
    (sub_lane(a.0, b.0), sub_lane(a.1, b.1), sub_lane(a.2, b.2), sub_lane(a.3, b.3))
}

#[inline(always)]
pub unsafe fn sub_const_8x(a: EF8, k: u64) -> EF8 {
    let q = _mm512_set1_epi64(4294967197);
    let k_vec = _mm512_set1_epi64(k as i64);
    let a0_plus_q = _mm512_add_epi64(a.0, q);
    let diff = _mm512_sub_epi64(a0_plus_q, k_vec);
    let m = _mm512_cmpge_epu64_mask(diff, q);
    let new_a0 = _mm512_mask_sub_epi64(diff, m, diff, q);
    (new_a0, a.1, a.2, a.3)
}

#[inline(always)]
pub unsafe fn ext_scalar_mul_8_8x(a: EF8) -> EF8 {
    let a2 = extadd_8x(a, a); let a4 = extadd_8x(a2, a2); extadd_8x(a4, a4)
}

#[inline(always)]
pub unsafe fn is_all_zeros_8x(a: EF8) -> bool {
    let z = _mm512_setzero_si512();
    let m0 = _mm512_cmpneq_epi64_mask(a.0, z); let m1 = _mm512_cmpneq_epi64_mask(a.1, z);
    let m2 = _mm512_cmpneq_epi64_mask(a.2, z); let m3 = _mm512_cmpneq_epi64_mask(a.3, z);
    (m0 | m1 | m2 | m3) == 0
}

#[inline(always)]
pub unsafe fn sum_soa_to_aos(soa: EF8) -> [u32; 4] {
    let sum_lane = |v: __m512i| -> u32 {
        let mut arr = [0u64; 8];
        _mm512_storeu_si512(arr.as_mut_ptr() as *mut _, v);
        let mut s = 0u64;
        for &val in arr.iter() { s += val; }
        (s % 4294967197) as u32
    };
    [sum_lane(soa.0), sum_lane(soa.1), sum_lane(soa.2), sum_lane(soa.3)]
}

#[inline(always)]
pub unsafe fn extmul_8x_local(a: EF8, b: EF8) -> EF8 {
    let mask32 = _mm512_set1_epi64(0xFFFFFFFF);
    let c99 = _mm512_set1_epi64(99);
    let red_simd = |x: __m512i| -> __m512i {
        let hi = _mm512_srli_epi64(x, 32); let lo = _mm512_and_si512(x, mask32);
        _mm512_add_epi64(lo, _mm512_mullo_epi64(hi, c99))
    };
    let p00 = red_simd(_mm512_mul_epu32(a.0, b.0)); let p01 = red_simd(_mm512_mul_epu32(a.0, b.1));
    let p02 = red_simd(_mm512_mul_epu32(a.0, b.2)); let p03 = red_simd(_mm512_mul_epu32(a.0, b.3));
    let p10 = red_simd(_mm512_mul_epu32(a.1, b.0)); let p11 = red_simd(_mm512_mul_epu32(a.1, b.1));
    let p12 = red_simd(_mm512_mul_epu32(a.1, b.2)); let p13 = red_simd(_mm512_mul_epu32(a.1, b.3));
    let p20 = red_simd(_mm512_mul_epu32(a.2, b.0)); let p21 = red_simd(_mm512_mul_epu32(a.2, b.1));
    let p22 = red_simd(_mm512_mul_epu32(a.2, b.2)); let p23 = red_simd(_mm512_mul_epu32(a.2, b.3));
    let p30 = red_simd(_mm512_mul_epu32(a.3, b.0)); let p31 = red_simd(_mm512_mul_epu32(a.3, b.1));
    let p32 = red_simd(_mm512_mul_epu32(a.3, b.2)); let p33 = red_simd(_mm512_mul_epu32(a.3, b.3));

    let t0 = p00; let t1 = _mm512_add_epi64(p01, p10); let t2 = _mm512_add_epi64(_mm512_add_epi64(p02, p11), p20);
    let t3 = _mm512_add_epi64(_mm512_add_epi64(_mm512_add_epi64(p03, p12), p21), p30);
    let t4 = _mm512_add_epi64(_mm512_add_epi64(p13, p22), p31); let t5 = _mm512_add_epi64(p23, p32); let t6 = p33;

    let final_red = |mut x: __m512i| -> __m512i {
        x = red_simd(x); x = red_simd(x);
        let q_vec = _mm512_set1_epi64(4294967197); let ge_q = _mm512_cmpge_epu64_mask(x, q_vec);
        _mm512_mask_sub_epi64(x, ge_q, x, q_vec)
    };

    (final_red(_mm512_add_epi64(t0, _mm512_slli_epi64(t4, 1))), final_red(_mm512_add_epi64(t1, _mm512_slli_epi64(t5, 1))),
     final_red(_mm512_add_epi64(t2, _mm512_slli_epi64(t6, 1))), final_red(t3))
}

#[inline(always)]
pub unsafe fn ext_fma_scalar_const_8x(z2: EF8, z: EF8, c_mod: u32, d: u32) -> EF8 {
    let c_vec = _mm512_set1_epi64(c_mod as i64); let d_vec = _mm512_set1_epi64(d as i64);
    let mask32 = _mm512_set1_epi64(0xFFFFFFFF); let c99 = _mm512_set1_epi64(99);
    let red_simd = |x: __m512i| -> __m512i {
        let hi = _mm512_srli_epi64(x, 32); let lo = _mm512_and_si512(x, mask32);
        _mm512_add_epi64(lo, _mm512_mullo_epi64(hi, c99))
    };
    let fma_lane = |z2_lane: __m512i, z_lane: __m512i, is_lane0: bool| -> __m512i {
        let p = red_simd(_mm512_mul_epu32(z_lane, c_vec)); 
        let mut sum = _mm512_add_epi64(z2_lane, p);
        if is_lane0 { sum = _mm512_add_epi64(sum, d_vec); }
        let mut res = red_simd(sum); res = red_simd(res);
        let q_vec = _mm512_set1_epi64(4294967197); let ge_q = _mm512_cmpge_epu64_mask(res, q_vec);
        _mm512_mask_sub_epi64(res, ge_q, res, q_vec)
    };
    (fma_lane(z2.0, z.0, true), fma_lane(z2.1, z.1, false), fma_lane(z2.2, z.2, false), fma_lane(z2.3, z.3, false))
}

#[inline(always)]
pub unsafe fn fold_8x(t0: EF8, t1: EF8, r_soa: EF8) -> EF8 {
    extadd_8x(t0, extmul_8x_local(extsub_8x(t1, t0), r_soa))
}