#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::macros::*;
// RNS: [257, 3329, 7681, 7937, 9473, 10753]

// @reduce 2^32-99
#[target_feature(enable = "avx2")]
pub unsafe fn reduce32(even_64: __m256i, odd_64: __m256i) -> __m256i {
    
    let reduce_lane = |x: __m256i| -> __m256i {
        let c99 = _mm256_set1_epi64x(99);
        let p_val = _mm256_set1_epi64x(4294967197); // 2^32 - 99
        let mask_low = _mm256_set1_epi64x(0x00000000FFFFFFFF);
        let zero = _mm256_setzero_si256();

        // Round 1
        let h1 = _mm256_srai_epi64::<32>(x); 
        let l1 = _mm256_and_si256(x, mask_low);
        let sum1 = _mm256_add_epi64(_mm256_mul_epi32(h1, c99), l1);

        // Round 2 (Convergence)
        let h2 = _mm256_srai_epi64::<32>(sum1);
        let l2 = _mm256_and_si256(sum1, mask_low);
        let sum2 = _mm256_add_epi64(_mm256_mul_epi32(h2, c99), l2);

        // Correction: if sum2 < 0, add P
        let is_neg = _mm256_cmpgt_epi64(zero, sum2);
        let res = _mm256_add_epi64(sum2, _mm256_and_si256(is_neg, p_val));

        _mm256_and_si256(res, mask_low)
    };

    let res_even = reduce_lane(even_64);
    let res_odd  = reduce_lane(odd_64);

    let odd_shifted = _mm256_slli_epi64::<32>(res_odd);
    _mm256_or_si256(res_even, odd_shifted)
}

// @barrett_fake - tested
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(9));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(7681));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_10753(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(6));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(10753));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_11777(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(6));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(11777));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_12289(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(5));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(12289));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_13313(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(5));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(13313));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_15361(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(4));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(15361));
    let f = _mm256_sub_epi16(x, e);
    f
}

// """ Sage
// def to_int16(x):
//     u = Integer(x) & 0xFFFF 
//     if u >= 0x8000:
//         return u - 0x10000
//     return u

// f = [769, 3329, 7681, 7937, 9473, 10753]
// for i in f:
//     for j in range(10000):
//         a = ZZ.random_element(-i, i)
//         r = ZZ.random_element(-(2^15), 2^15-1)
//         m = ZZ.random_element(-(2^15), 2^15-1)
//         m1 = round((m*(2^16))/(i))
//         t = (r*m1) >> 16
//         d = to_int16(r*m)
//         dp = to_int16(a+d)
//         ds = to_int16(a-d)
//         ansp = to_int16(dp - to_int16(t*i))
//         anss = to_int16(ds + to_int16(t*i))
//         if ansp%i != (a + r*m)%i or ansp - (a + (r*m)%i) < -3*i or ansp - (a + (r*m)%i) > 3*i:
//             print(ansp, a + (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
//         if anss%i != (a - r*m)%i or anss - (a - (r*m)%i) < -3*i or anss - (a - (r*m)%i) > 3*i:
//             print(anss, a - (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
// """

// @barrett_butterfly - tested
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(7681));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 7681);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 7681);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(10753));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 10753);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 10753);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_11777(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(11777));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 11777);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 11777);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_12289(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(12289));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 12289);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 12289);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_13313(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(13313));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 13313);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 13313);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_15361(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_mullo_epi16(t, _mm256_set1_epi16(15361));
    let ansp = range_narrow(_mm256_sub_epi16(dp, e), 15361);
    let anss = range_narrow(_mm256_add_epi16(ds, e), 15361);
    [ansp, anss]
}


// @barrett_fake_32 
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(5585133);
    let m2 = _mm256_set1_epi32(7681);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(5585133));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(7681));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(5585133));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(7681));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_10753_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(1290167);
    let m2 = _mm256_set1_epi32(10753);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(1290167));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(10753));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(1290167));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(10753));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_11777_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(559168);
    let m2 = _mm256_set1_epi32(11777);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(559168));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(11777));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(559168));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(11777));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_12289_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(541132);
    let m2 = _mm256_set1_epi32(12289);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(541132));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(12289));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(541132));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(12289));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_13313_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(453390);
    let m2 = _mm256_set1_epi32(13313);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(453390));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(13313));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(453390));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(13313));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_15361_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(399420);
    let m2 = _mm256_set1_epi32(15361);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(399420));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(15361));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(399420));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(15361));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
    ans
}

#[inline]
pub unsafe fn range_narrow(mut input: __m256i, q: i16) -> __m256i {
    let q_vec = _mm256_set1_epi16(q);
    let half_q = _mm256_set1_epi16(q / 2);
    let neg_half_q = _mm256_set1_epi16(-(q / 2));

    let up_mask1 = _mm256_cmpgt_epi16(neg_half_q, input);
    let down_mask1 = _mm256_cmpgt_epi16(input, half_q);

    input = _mm256_add_epi16(input, _mm256_and_si256(up_mask1, q_vec));
    input = _mm256_sub_epi16(input, _mm256_and_si256(down_mask1, q_vec));

    let up_mask2 = _mm256_cmpgt_epi16(neg_half_q, input);
    let down_mask2 = _mm256_cmpgt_epi16(input, half_q);

    input = _mm256_add_epi16(input, _mm256_and_si256(up_mask2, q_vec));
    input = _mm256_sub_epi16(input, _mm256_and_si256(down_mask2, q_vec));

    input
}

// @barrett_mul - tested
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7681(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 7681);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(7681));
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 7681)
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_10753(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 10753);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(10753));
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 10753)
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_11777(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 11777);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(11777));
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 11777)
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_12289(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 12289);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(12289));
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 12289)
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_13313(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 13313);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(13313));
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 13313)
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_15361(bb: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let b = range_narrow(bb, 15361);
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(15361));       
    let g = _mm256_sub_epi16(d, f);
    range_narrow(g, 15361)
}

// @barrett_ibutterfly - tested
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_7681(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_7681(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_10753(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_10753(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_11777(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_11777(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_11777(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_12289(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_12289(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_12289(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_13313(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_13313(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_13313(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_15361(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i;2]{
    let ans0 = barrett_fake_15361(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul_15361(tmp, c, cr);
    [ans0, ans1]
}

// @montproduct

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_7681(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(7681);
    let mu_vec = _mm256_set1_epi16(-7679);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 7681)
}

// --- Field: 10753 ---
#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_10753(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(10753);
    let mu_vec = _mm256_set1_epi16(-10751);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 10753)
}

// --- Field: 11777 ---
#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_11777(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(11777);
    let mu_vec = _mm256_set1_epi16(-11775);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 11777)
}

// --- Field: 12289 ---
#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_12289(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(12289);
    let mu_vec = _mm256_set1_epi16(-12287);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 12289)
}

// --- Field: 13313 ---
#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_13313(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(13313);
    let mu_vec = _mm256_set1_epi16(-13311);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 13313)
}

// --- Field: 15361 ---
#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_15361(x: __m256i, y: __m256i) -> __m256i {
    let f_vec = _mm256_set1_epi16(15361);
    let mu_vec = _mm256_set1_epi16(-15359);
    
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    
    // k = (x * y * mu) mod 2^16
    let k = _mm256_mullo_epi16(lo, mu_vec);
    
    // res = (x * y - k * f) / 2^16
    let kf_hi = _mm256_mulhi_epi16(k, f_vec);
    let res = _mm256_sub_epi16(hi, kf_hi);
    range_narrow(res, 15361)
}




#[cfg(test)]
mod tests {
    use super::*; 
    use core::arch::x86_64::*;

    unsafe fn dump_m256i(v: __m256i) -> [i16; 16] {
        let mut arr = [0i16; 16];
        
        _mm256_storeu_si256(arr.as_mut_ptr() as *mut __m256i, v);
        arr
    }

    macro_rules! test_barrett {
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;
                    let qq = $q as i16;
                    let mut inputs = [0i16; 16];
                    
                    for start_val in (-32768..32768).step_by(16) {
                        for i in 0..16 {
                            inputs[i] = (start_val + i as i32) as i16;
                        }

                        let input_vec = _mm256_loadu_si256(inputs.as_ptr() as *const __m256i);
                        let output_vec = super::$func_name(input_vec);
                        let outputs = dump_m256i(output_vec);

                        for i in 0..16 {
                            let x = inputs[i];
                            let y = outputs[i];
                            
                            let diff = (x as i32) - (y as i32);
                            if diff % q != 0 {
                                println!("\n=== EXHAUSTIVE TEST FAILED: {} (q={}) ===", stringify!($func_name), q);
                                println!("Input x:   {}", x);
                                println!("Output y:  {}", y);
                                println!("Diff:      {}", diff);
                                println!("Modulo q:  {}", diff % q);
                                
                                panic!(
                                    "Failed for {}: x={}, y={}", 
                                    stringify!($func_name), x, y
                                );
                            }

                            let bound = qq+3000;

                            if y < -bound || y > bound {
                                panic!(
                                    "\n=== RANGE ERROR: {} (q={}) ===\nInput x: {}\nOutput y: {}\nOutput is outside the range (-q, q)!", 
                                    stringify!($func_name), q, x, y
                                );
                            }
                        }
                    }
                    println!("SUCCESS: {} passed exhaustive test for all 65536 values.", stringify!($func_name));
                }
            }
        };
    }

    // test_barrett!(barrett_fake_7681, 7681);
    // test_barrett!(barrett_fake_10753, 10753);
    // test_barrett!(barrett_fake_11777, 11777);
    // test_barrett!(barrett_fake_12289, 12289);
    // test_barrett!(barrett_fake_13313, 13313);
    // test_barrett!(barrett_fake_15361, 15361);

    macro_rules! test_butterfly {
        ($func_name:ident, $q:expr, $c_list:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;
                    let qq = $q as i16;
                    let c_test_cases: &[i16] = &$c_list;

                    for &c_i16 in c_test_cases {
                        let c = c_i16 as i32;
                        let cr_val = ((c_i16 as f64 * 65536.0) / q as f64).round() as i16;
                        
                        let c_vec = _mm256_set1_epi16(c_i16);
                        let cr_vec = _mm256_set1_epi16(cr_val);

                        for a_base in (-q-300..q+300).step_by(256) {
                            for b_base in (-q-300..q+300).step_by(256) {
                                let mut a_vals = [0i16; 16];
                                let mut b_vals = [0i16; 16];
                                
                                for i in 0..16 {
                                    a_vals[i] = (a_base as i16).wrapping_add(i as i16);
                                    b_vals[i] = (b_base as i16).wrapping_add((i * 7) as i16);
                                }

                                let a_vec = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
                                let b_vec = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

                                let [ansp_vec, anss_vec] = super::$func_name(a_vec, b_vec, c_vec, cr_vec);
                                
                                let ansp_res = dump_m256i(ansp_vec);
                                let anss_res = dump_m256i(anss_vec);

                                for i in 0..16 {
                                    let a_val = a_vals[i] as i32;
                                    let b_val = b_vals[i] as i32;

                                    let expected_p = (a_val + b_val * c).rem_euclid(q);
                                    let expected_s = (a_val - b_val * c).rem_euclid(q);

                                    let res_p = (ansp_res[i] as i32).rem_euclid(q);
                                    let res_s = (anss_res[i] as i32).rem_euclid(q);

                                    assert_eq!(res_p, expected_p, 
                                        "\nFAILED P: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a_val, b_val, ansp_res[i], expected_p);

                                    assert_eq!(res_s, expected_s, 
                                        "\nFAILED S: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a_val, b_val, anss_res[i], expected_s);

                                    assert!(ansp_res[i] > (-qq) && (ansp_res[i] < qq),
                                        "\nFAILED range P: q={}, c={}\na={}, b={}\nGot: {}", 
                                        q, c, a_val, b_val, ansp_res[i]);

                                    assert!(anss_res[i] > (-qq) && (anss_res[i] < qq),
                                        "\nFAILED range S: q={}, c={}\na={}, b={}\nGot: {}", 
                                        q, c, a_val, b_val, anss_res[i]);
                                }
                            }
                        }
                    }
                    println!("SUCCESS: {} passed for all c cases and i16 range.", stringify!($func_name));
                }
            }
        };
    }

    // test_butterfly!(barrett_butterfly_7681, 7681, [1925, 1213, 1213, 1925, 583, 527, 1728, 849, 527, 583, 849, 1728, 2381, 2784, 2446, 1366, 97, -2138, 2132, 2648, 2784, 2381, 1366, 2446, 2138, -97, 2648, 2132, 2381, -2381, 2446, -2446, 97, -97, 2132, -2132, 2784, -2784, 1366, -1366, 2138, -2138, 2648, -2648, 878, 2273, 330, -2645, -2753, -3654, 365, -1846, 1286, 3092, 675, -2268, -1112, -1794, 2399, -3000, 1286, -1286, 3092, -3092, 675, -675, 2268, -2268, 1794, -1794, 1112, -1112, 2399, -2399, 3000, -3000, 878, -878, 2273, -2273, 330, -330, 2645, -2645, 3654, -3654, 2753, -2753, 365, -365, 1846, -1846, -799, -695, -1381, 1875, 2724, 1908, 1382, -2423, -2469, -3380, 693, -1714, -3080, 3477, -732, 3074, -1740, -2774, -584, 1655, 2941, 2508, -3449, 528, -1080, 2516, -3411, -2551, -202, 243, 2881, -766, 1740, -1740, 2774, -2774, 584, -584, 1655, -1655, 2941, -2941, 2508, -2508, 528, -528, 3449, -3449, 2551, -2551, 3411, -3411, 2516, -2516, 1080, -1080, 202, -202, 243, -243, 766, -766, 2881, -2881, 799, -799, 695, -695, 1381, -1381, 1875, -1875, 2724, -2724, 1908, -1908, 2423, -2423, 1382, -1382, 1714, -1714, 693, -693, 3380, -3380, 2469, -2469, 3080, -3080, 3477, -3477, 3074, -3074, 732, -732, -1587, 198, 2900, -2063, -880, -3188, 219, -3501, 405, -2897, 319, 3837, 1633, 1800, 1996, -869, -1155, 2264, -3566, 3073, -1886, 2573, -1220, -2563, -257, 1478, 3141, 3180, -1402, 3789, 3125, -2819, 1125, -3780, -417, 2593, -2990, -707, 1438, -2681, 550, -1848, 1228, 1097, -1591, 2028, 1952, 2044, -2880, -3532, -1682, 1415, 2562, -3078, 648, -3099, 1003, -1853, 2844, 3041, 993, -2722, 1044, -1408, 3780, -3780, 1125, -1125, 2593, -2593, 417, -417, 707, -707, 2990, -2990, 1438, -1438, 2681, -2681, 550, -550, 1848, -1848, 1228, -1228, 1097, -1097, 1952, -1952, 2044, -2044, 2028, -2028, 1591, -1591, 3099, -3099, 648, -648, 3078, -3078, 2562, -2562, 1682, -1682, 1415, -1415, 3532, -3532, 2880, -2880, 1853, -1853, 1003, -1003, 2844, -2844, 3041, -3041, 1408, -1408, 1044, -1044, 993, -993, 2722, -2722, 198, -198, 1587, -1587, 2063, -2063, 2900, -2900, 3188, -3188, 880, -880, 219, -219, 3501, -3501, 405, -405, 2897, -2897, 319, -319, 3837, -3837, 869, -869, 1996, -1996, 1633, -1633, 1800, -1800, 2563, -2563, 1220, -1220, 1886, -1886, 2573, -2573, 3073, -3073, 3566, -3566, 2264, -2264, 1155, -1155, 1478, -1478, 257, -257, 3141, -3141, 3180, -3180, 3125, -3125, 2819, -2819, 3789, -3789, 1402, -1402, 763, -413, -1704, -3799, -2689, 2583, -2668, -669, 2358, -3445, 2922, 321, 1656, -2799, 3694, 185, -1950, 1129, -398, 2259, -3546, -1604, -2359, 62, 1607, 1667, -1968, 1683, 201, -3626, -1979, 2875, 3081, 94, -1193, -3394, 1131, 1035, 3452, -2996, -2173, 542, -3120, -1266, -1437, 702, 1065, 506, -1230, -2012, -1876, -2002, 2197, 2757, 3006, 346, 3139, -3586, -2169, -2372, 296, 2838, 1959, -1406, -3250, -3239, -1897, 3765, -2457, -1189, -1771, 113, 335, 3483, -738, 329, 218, -118, 3280, 2805, -2840, 1211, 3832, -1872, -639, -3376, -674, 1115, 3751, -621, 535, -2811, 3016, 2760, -2252, 1036, 2110, 2481, 1499, -1657, 2395, -1170, -1775, -1717, 3193, -2433, 1885, 1725, -2717, -2546, 572, 536, -217, -3265, -2951, -2067, 3615, 1393, 856, 111, 2050, 793, -1994, 1784, -1459, -3086, -2671, 3137] );
    // test_butterfly!(barrett_butterfly_10753, 10753, [4489, 67, 321, 321, 67, 3422, 1656, 4679, 3461, 1656, 3422, 3461, 4679, 1560, 2640, 2637, -1154, 4631, -4832, 3010, 2047, 2640, 1560, 1154, -2637, 4832, -4631, 2047, 3010, 1560, -1560, 2637, -2637, 4631, -4631, 3010, -3010, 2640, -2640, 1154, -1154, 4832, -4832, 2047, -2047, -4611, -746, -3783, -2900, -4351, 4191, 1219, 1186, -597, -2436, -1917, -3013, 644, -1641, 2417, 136, 597, -597, 2436, -2436, 3013, -3013, 1917, -1917, 1641, -1641, 644, -644, 2417, -2417, 136, -136, 4611, -4611, 746, -746, 2900, -2900, 3783, -3783, 4191, -4191, 4351, -4351, 1219, -1219, 1186, -1186, -3907, -380, -3954, 3697, 5147, 3314, 753, 3775, -3982, -3712, -1385, 2031, 2603, -3644, 3171, 2353, -1047, -922, 5122, -2744, 94, 2599, 4455, 2085, -559, 3902, -5194, -3362, 946, -841, 1136, 2582, 922, -922, 1047, -1047, 2744, -2744, 5122, -5122, 2085, -2085, 4455, -4455, 2599, -2599, 94, -94, 3362, -3362, 5194, -5194, 559, -559, 3902, -3902, 946, -946, 841, -841, 1136, -1136, 2582, -2582, 380, -380, 3907, -3907, 3697, -3697, 3954, -3954, 3775, -3775, 753, -753, 3314, -3314, 5147, -5147, 1385, -1385, 2031, -2031, 3712, -3712, 3982, -3982, 2603, -2603, 3644, -3644, 3171, -3171, 2353, -2353, -685, 393, -2883, 4825, -5125, -5295, 721, -84, -2004, 4305, -1896, 5232, -2726, -100, 159, -4053, 267, -4980, 3617, -317, -671, -1279, 331, 1945, 3719, 4818, 1854, 216, 118, -2805, -5134, 2847, 2135, -3092, 3256, -2857, 4683, -128, -2177, -1924, 4193, -4627, -1353, -1828, -5178, 3944, -4577, 2830, -5263, -1266, 1202, 2228, -1289, 1207, -339, 5155, 1943, 1444, 1145, 29, 5012, -3592, -2461, 4098, 2857, -2857, 3256, -3256, 3092, -3092, 2135, -2135, 1924, -1924, 2177, -2177, 4683, -4683, 128, -128, 4577, -4577, 2830, -2830, 5178, -5178, 3944, -3944, 1353, -1353, 1828, -1828, 4193, -4193, 4627, -4627, 5155, -5155, 339, -339, 1207, -1207, 1289, -1289, 1266, -1266, 5263, -5263, 1202, -1202, 2228, -2228, 1943, -1943, 1444, -1444, 29, -29, 1145, -1145, 5012, -5012, 3592, -3592, 2461, -2461, 4098, -4098, 4825, -4825, 2883, -2883, 393, -393, 685, -685, 721, -721, 84, -84, 5295, -5295, 5125, -5125, 159, -159, 4053, -4053, 2726, -2726, 100, -100, 1896, -1896, 5232, -5232, 2004, -2004, 4305, -4305, 1279, -1279, 671, -671, 331, -331, 1945, -1945, 317, -317, 3617, -3617, 4980, -4980, 267, -267, 3719, -3719, 4818, -4818, 216, -216, 1854, -1854, 118, -118, 2805, -2805, 5134, -5134, 2847, -2847, -1102, -498, 1437, -1107, 3259, -5182, -3293, -3098, -2129, 2336, 4783, 2854, 5096, -4313, -2664, -1360, -1036, -5308, 4894, -787, 4847, -4864, -2159, 3298, 670, 3210, 1878, 10, -2351, 4946, -1961, 3778, 1196, 3097, -3192, -4861, -4524, -4181, -2024, 549, -4711, 3472, -3800, -3942, 3223, -5262, -881, 2295, -2545, 4819, 283, 1533, -940, -4484, -656, 1538, -3992, -5163, -1825, 1361, 2343, 1293, 607, -4314, -1907, -1135, -774, 1267, 1317, -2137, -3390, -2215, 3661, -3645, -2032, -3104, 2076, 3687, 290, -697, -2758, 3959, 3572, -1985, 1082, 3258, 3226, 2777, -2266, 264, -1280, 3818, -3689, -301, 156, -1339, -4209, -1180, 2425, -3789, 2966, -2160, -4931, 5168, 2037, -4043, -2056, -3310, -2670, -3965, 3911, 3170, -1445, 2546, 38, 1466, -1590, 2482, -4999, -1000, -3429, -5238, 3903, -3930, 815, 2515, 3543, -840] );
    // test_butterfly!(barrett_butterfly_11777, 11777, [5322, 4976, 4201, 4201, 4976, 2497, -3410, 4578, -337, 3410, -2497, 337, -4578, 1224, -4782, 1447, -293, 1915, 2380, 4525, 5692, 293, -1447, 4782, -1224, 5692, 4525, 2380, 1915, 1224, -1224, 1447, -1447, 1915, -1915, 4525, -4525, 293, -293, 4782, -4782, 5692, -5692, 2380, -2380, -4633, 4212, -4148, 5558, 5573, -5020, 503, -3587, 1265, -4114, 2838, -5722, 3114, -2469, 2333, -3268, 5722, -5722, 2838, -2838, 4114, -4114, 1265, -1265, 2333, -2333, 3268, -3268, 3114, -3114, 2469, -2469, 5558, -5558, 4148, -4148, 4212, -4212, 4633, -4633, 503, -503, 3587, -3587, 5573, -5573, 5020, -5020, 1345, -2326, 3384, -2615, -5739, -5197, -2020, 1961, -1548, -5444, -2244, -690, 3488, -2584, 2500, -3010, -309, 4282, -5203, -2639, 5537, -1860, -1362, 5709, 3238, -2885, -403, 1352, -4055, -5246, -5513, -3679, 5709, -5709, 1362, -1362, 5537, -5537, 1860, -1860, 2639, -2639, 5203, -5203, 309, -309, 4282, -4282, 4055, -4055, 5246, -5246, 3679, -3679, 5513, -5513, 3238, -3238, 2885, -2885, 1352, -1352, 403, -403, 2020, -2020, 1961, -1961, 5197, -5197, 5739, -5739, 3384, -3384, 2615, -2615, 2326, -2326, 1345, -1345, 3488, -3488, 2584, -2584, 3010, -3010, 2500, -2500, 1548, -1548, 5444, -5444, 690, -690, 2244, -2244, -1386, 3890, 4586, -4748, 1596, -2695, -3683, -3998, -1120, -1478, 2599, -5680, -3452, -576, 5491, -4365, 2873, 3560, -1952, 1230, 1534, 2487, 1688, 2315, 4697, -5137, 5622, 5073, -50, 4771, -1483, -1936, -1329, -5039, 5570, -831, -4547, 2601, -2255, -347, 4302, -756, -4993, -3834, -1203, 4322, -1470, 3412, 4779, -4482, 3206, 2541, -2926, -2978, -3062, 3404, 2802, 2562, 1216, -5798, 3673, 2114, 1056, 2403, 1470, -1470, 3412, -3412, 4322, -4322, 1203, -1203, 4302, -4302, 756, -756, 3834, -3834, 4993, -4993, 347, -347, 2255, -2255, 2601, -2601, 4547, -4547, 5039, -5039, 1329, -1329, 5570, -5570, 831, -831, 2562, -2562, 2802, -2802, 5798, -5798, 1216, -1216, 2403, -2403, 1056, -1056, 2114, -2114, 3673, -3673, 4779, -4779, 4482, -4482, 2541, -2541, 3206, -3206, 3062, -3062, 3404, -3404, 2978, -2978, 2926, -2926, 576, -576, 3452, -3452, 5491, -5491, 4365, -4365, 5680, -5680, 2599, -2599, 1478, -1478, 1120, -1120, 1596, -1596, 2695, -2695, 3998, -3998, 3683, -3683, 4748, -4748, 4586, -4586, 1386, -1386, 3890, -3890, 4697, -4697, 5137, -5137, 5073, -5073, 5622, -5622, 1936, -1936, 1483, -1483, 50, -50, 4771, -4771, 3560, -3560, 2873, -2873, 1230, -1230, 1952, -1952, 2315, -2315, 1688, -1688, 2487, -2487, 1534, -1534, 4905, 5199, -3845, 5336, 295, -3649, -2710, 4205, 2550, 3996, 4520, -4971, -371, -4074, -2893, 4007, -1148, 2607, 603, 5822, 3002, 4745, 1765, 4716, -5169, 1654, 1819, 24, 599, 3689, -3879, 1043, -2105, -2883, -4727, -1422, 5857, -2765, 3104, -3643, 3225, 4361, 4451, -4675, 2643, -4308, -2468, 3341, -3131, 1273, 1115, 1578, 5029, 4783, -1109, 1821, 3875, 1223, -3051, 3061, -4819, -3588, -1372, 44, 3435, 3166, 3610, -4133, 3539, 3135, 4765, -3449, 2265, 5318, 51, 551, 2745, 5410, 2200, -2062, 5713, 3628, -1213, 1790, 2172, -5630, -3414, -2603, 5495, -2099, 3074, -1575, 743, 2826, -438, 810, -5518, 5042, 5381, 3982, -1200, 3266, 644, 261, 803, -1485, -3325, -5181, 260, 5811, -1710, 3001, -4071, -3818, -856, -2067, 5824, -1736, -2973, 5795, -3326, 141, 5004, -3491, 1233, -2237, 409, 2047] );
    // test_butterfly!(barrett_butterfly_12289, 12289, [1479, 4043, 5146, 5146, 4043, 5736, -722, 4134, -1305, 722, -5736, 1305, -4134, 1646, 3621, 1212, 2545, 5860, -3542, 3195, 3504, 2545, 1212, 3621, 1646, 3504, 3195, 3542, -5860, 1646, -1646, 1212, -1212, 5860, -5860, 3195, -3195, 2545, -2545, 3621, -3621, 3504, -3504, 3542, -3542, -3006, 2744, 2975, 563, 2639, 4821, -949, -2625, 5023, 5828, -4591, 5728, -3328, -5777, 1351, 4978, 4591, -4591, 5728, -5728, 5023, -5023, 5828, -5828, 4978, -4978, 1351, -1351, 3328, -3328, 5777, -5777, 2975, -2975, 563, -563, 3006, -3006, 2744, -2744, 949, -949, 2625, -2625, 4821, -4821, 2639, -2639, 4805, -3553, 2294, -1062, 3712, -3135, -2747, -4846, -955, -790, -1170, 2319, 5086, -1326, -3201, 3014, 2963, -4896, -3051, -2366, -4320, 1000, 3091, -81, -1177, 4255, 1635, -2768, -4611, 726, -1853, -140, 1000, -1000, 4320, -4320, 3091, -3091, 81, -81, 2963, -2963, 4896, -4896, 3051, -3051, 2366, -2366, 1853, -1853, 140, -140, 4611, -4611, 726, -726, 4255, -4255, 1177, -1177, 2768, -2768, 1635, -1635, 3712, -3712, 3135, -3135, 2747, -2747, 4846, -4846, 3553, -3553, 4805, -4805, 2294, -2294, 1062, -1062, 1326, -1326, 5086, -5086, 3014, -3014, 3201, -3201, 1170, -1170, 2319, -2319, 955, -955, 790, -790, -3637, -3459, -145, -5542, -5911, 4890, 2731, 3932, 2548, -4231, 3382, -355, -1759, -3707, 5179, -3694, -2426, -334, 1428, -1696, -4632, -5755, 1260, -4388, -2881, -3284, -5092, -2089, 2013, 3289, -729, -3241, 4278, -1673, -5331, 4989, -4177, -3584, -2525, 1381, -2468, -339, -544, -5791, -2842, 480, 1022, 9, -953, -3748, -5767, -827, -118, -2476, 2197, 5067, -4354, -130, 5374, -2837, -4452, -2396, -3949, 3296, 544, -544, 5791, -5791, 339, -339, 2468, -2468, 2842, -2842, 480, -480, 9, -9, 1022, -1022, 4278, -4278, 1673, -1673, 4989, -4989, 5331, -5331, 3584, -3584, 4177, -4177, 1381, -1381, 2525, -2525, 2396, -2396, 4452, -4452, 3296, -3296, 3949, -3949, 130, -130, 4354, -4354, 5374, -5374, 2837, -2837, 5767, -5767, 827, -827, 3748, -3748, 953, -953, 5067, -5067, 2197, -2197, 118, -118, 2476, -2476, 2548, -2548, 4231, -4231, 355, -355, 3382, -3382, 3707, -3707, 1759, -1759, 3694, -3694, 5179, -5179, 5542, -5542, 145, -145, 3637, -3637, 3459, -3459, 5911, -5911, 4890, -4890, 3932, -3932, 2731, -2731, 2089, -2089, 5092, -5092, 2881, -2881, 3284, -3284, 729, -729, 3241, -3241, 3289, -3289, 2013, -2013, 5755, -5755, 4632, -4632, 1260, -1260, 4388, -4388, 334, -334, 2426, -2426, 1696, -1696, 1428, -1428, -2143, -1065, -4645, -404, -5277, -1168, 1207, -3248, -1378, -1912, -435, -4337, -4096, -493, 2381, 5444, 2187, -2566, 6039, 2422, 2437, -3646, 6022, 2987, 4976, 1607, -3780, -875, 4284, 5088, -1002, 5011, -4861, -354, 5698, 2912, 5012, -2481, -1045, 2859, 773, -390, 3833, 3778, -2401, -442, 1067, 5101, -242, -1537, -4143, -4714, -2678, 3704, 5019, -545, 1017, -4885, -5084, 1632, 3066, -27, 1440, -3763, 4057, -3271, 3364, 1689, 2847, 4414, 4372, -2174, -5042, -2305, 4053, 2645, 4895, -1484, 5195, 2780, -3636, 4938, -2704, -5291, -1777, -1663, -4654, -1426, -2166, -3915, -113, -4919, 160, 3149, -3, -4437, -3247, 2686, 3978, -2969, -3510, -5332, -2865, -2370, -3186, 5407, -2126, 1630, -1153, 2884, 4048, -2249, 5559, 420, -1544, 2178, -476, -3531, 4905, 3985, 671, 3000, -243, -3016, 3400, -2399, 3136, -5191] );
    // test_butterfly!(barrett_butterfly_13313, 13313, [258, 1455, -2626, 2626, -1455, 5275, -6630, 3024, 6476, 6476, 3024, 6630, -5275, 1415, 4215, 5619, 4196, 4690, 5487, 1463, 4468, 5487, 4690, 4468, 1463, 4196, 5619, 4215, 1415, 1415, -1415, 5619, -5619, 4690, -4690, 1463, -1463, 5487, -5487, 4468, -4468, 4196, -4196, 4215, -4215, -3101, 1278, -4330, -1152, 6072, -4358, 5072, 3902, 1173, -3565, 4995, 2651, -5894, -2970, -5375, 2198, 5375, -5375, 2198, -2198, 5894, -5894, 2970, -2970, 1173, -1173, 3565, -3565, 4995, -4995, 2651, -2651, 6072, -6072, 4358, -4358, 5072, -5072, 3902, -3902, 4330, -4330, 1152, -1152, 3101, -3101, 1278, -1278, -772, 519, 3696, -4968, -6168, 6216, 4753, 1478, 1229, -2430, -4253, 5608, -2169, -456, 714, -2170, 3813, -1408, 3606, -1562, -2627, -1197, 2368, -1454, 272, -3611, 3630, 4630, 6105, -4156, 2878, -3004, 4156, -4156, 6105, -6105, 2878, -2878, 3004, -3004, 3611, -3611, 272, -272, 4630, -4630, 3630, -3630, 3813, -3813, 1408, -1408, 1562, -1562, 3606, -3606, 2627, -2627, 1197, -1197, 2368, -2368, 1454, -1454, 1229, -1229, 2430, -2430, 5608, -5608, 4253, -4253, 2169, -2169, 456, -456, 714, -714, 2170, -2170, 6216, -6216, 6168, -6168, 1478, -1478, 4753, -4753, 519, -519, 772, -772, 3696, -3696, 4968, -4968, 3336, 4657, 5365, -382, -4787, -3063, 2386, 3190, 5893, 2712, 743, 5312, 5303, 3065, 5675, -280, 3042, 631, 6194, -492, 275, -4385, -3248, -735, 3904, -4556, 4331, -894, 2935, 1611, -3048, -917, 3523, 3650, -460, 1137, -3152, -1123, -6488, 3534, -225, 4798, -5078, 5450, 694, 5983, -2018, -1437, 2401, 6253, 5336, 5449, 5039, -4612, 3718, 708, 3345, -2335, 2610, 5583, -1969, -2108, 2600, 5150, 5583, -5583, 2610, -2610, 3345, -3345, 2335, -2335, 1969, -1969, 2108, -2108, 5150, -5150, 2600, -2600, 5449, -5449, 5336, -5336, 2401, -2401, 6253, -6253, 3718, -3718, 708, -708, 5039, -5039, 4612, -4612, 3523, -3523, 3650, -3650, 1137, -1137, 460, -460, 3534, -3534, 6488, -6488, 3152, -3152, 1123, -1123, 4798, -4798, 225, -225, 5450, -5450, 5078, -5078, 694, -694, 5983, -5983, 1437, -1437, 2018, -2018, 3042, -3042, 631, -631, 492, -492, 6194, -6194, 3248, -3248, 735, -735, 4385, -4385, 275, -275, 4556, -4556, 3904, -3904, 894, -894, 4331, -4331, 2935, -2935, 1611, -1611, 917, -917, 3048, -3048, 743, -743, 5312, -5312, 2712, -2712, 5893, -5893, 5675, -5675, 280, -280, 5303, -5303, 3065, -3065, 5365, -5365, 382, -382, 4657, -4657, 3336, -3336, 4787, -4787, 3063, -3063, 3190, -3190, 2386, -2386, -2518, -2693, -4293, 2615, -604, 3924, 162, 1857, 2909, 4994, 2628, 939, 3077, 4914, 3867, -789, -407, -1498, 6413, 3742, -3532, -5972, -242, -4129, 5405, 3375, -1872, -3708, 3644, -5071, -3446, 2903, 5730, 597, 3290, 3212, -5240, 6006, -5402, -4149, -4200, 5247, -6036, -333, -4796, -741, -198, 2168, -1600, 97, -1775, 5308, -5781, 442, -2461, 4086, 2288, -4532, -4125, -790, 5691, -3848, -281, -5933, 4282, 223, 174, 4953, -4781, -4611, -6344, 747, -5857, -6576, 1615, -3967, -3798, -5278, -1195, 2111, 3870, -15, -4801, 549, 6259, -3949, -5421, 753, 3625, -3340, -2427, -455, -4407, -5401, -4681, 3785, -49, -671, -4456, -4730, 5358, 2188, -5528, -1733, 3489, -5122, 4242, -2770, 6473, 5909, -5924, -2603, -3196, 838, 3943, 5506, -4816, 4419, -534, -4642, 5121, 3231, 4225, 1616, 2885, -1198, 4080, 913] );
    // test_butterfly!(barrett_butterfly_15361, 15361, [3968, 1319, 4309, 4309, 1319, 179, 3261, 3666, -5686, 3261, 179, 5686, -3666, 3874, 6372, 4329, -110, 5407, -6841, 4341, -2201, 6372, 3874, 110, -4329, 2201, -4341, 6841, -5407, 3874, -3874, 4329, -4329, 5407, -5407, 4341, -4341, 6372, -6372, 110, -110, 2201, -2201, 6841, -6841, 2395, 5099, -5361, 2537, 7237, 6707, -6422, -1403, 186, 720, 442, -2702, -7467, -2313, 2572, 5992, 186, -186, 720, -720, 2702, -2702, 442, -442, 5992, -5992, 2572, -2572, 2313, -2313, 7467, -7467, 2395, -2395, 5099, -5099, 2537, -2537, 5361, -5361, 1403, -1403, 6422, -6422, 6707, -6707, 7237, -7237, -3519, 243, -2064, -2539, -100, 2586, -6349, 792, 2792, 3375, 3992, -3065, -5046, 7145, -7399, 4361, -5989, 885, -3937, -121, -4805, 3239, 1883, -6298, 1733, -5184, 2051, -2962, 6276, -2987, 1535, 7437, 5989, -5989, 885, -885, 3937, -3937, 121, -121, 6298, -6298, 1883, -1883, 4805, -4805, 3239, -3239, 1535, -1535, 7437, -7437, 6276, -6276, 2987, -2987, 2962, -2962, 2051, -2051, 5184, -5184, 1733, -1733, 3519, -3519, 243, -243, 2064, -2064, 2539, -2539, 6349, -6349, 792, -792, 2586, -2586, 100, -100, 4361, -4361, 7399, -7399, 7145, -7145, 5046, -5046, 3992, -3992, 3065, -3065, 2792, -2792, 3375, -3375, -1100, -2276, 6649, -6966, -7374, -2673, 2793, 7343, 6403, -10, -2993, -2171, -1790, -5938, 1888, 4584, -3702, 4420, -1860, 7200, -7592, -2135, -5002, -1524, 2776, 1331, 4435, -5626, 4907, -6772, 7527, 5352, -5757, -1969, 1102, -5149, -1316, -852, 2435, 11, -6713, -1210, 1554, 6511, -1668, -1967, -1536, 3469, 2802, 3028, -6163, 72, 2473, -2815, 5355, 4377, -2430, -4468, -5332, -5279, -1000, 4862, -7441, 2046, 1969, -1969, 5757, -5757, 1102, -1102, 5149, -5149, 852, -852, 1316, -1316, 11, -11, 2435, -2435, 3469, -3469, 1536, -1536, 1668, -1668, 1967, -1967, 1210, -1210, 6713, -6713, 1554, -1554, 6511, -6511, 1000, -1000, 4862, -4862, 7441, -7441, 2046, -2046, 2430, -2430, 4468, -4468, 5279, -5279, 5332, -5332, 4377, -4377, 5355, -5355, 2473, -2473, 2815, -2815, 72, -72, 6163, -6163, 2802, -2802, 3028, -3028, 2276, -2276, 1100, -1100, 6649, -6649, 6966, -6966, 2673, -2673, 7374, -7374, 7343, -7343, 2793, -2793, 5938, -5938, 1790, -1790, 1888, -1888, 4584, -4584, 2993, -2993, 2171, -2171, 10, -10, 6403, -6403, 7527, -7527, 5352, -5352, 6772, -6772, 4907, -4907, 4435, -4435, 5626, -5626, 1331, -1331, 2776, -2776, 7592, -7592, 2135, -2135, 1524, -1524, 5002, -5002, 3702, -3702, 4420, -4420, 1860, -1860, 7200, -7200, 6440, 6784, 273, -7367, -2353, 2784, -685, 817, -1794, 6449, 3763, -692, -2296, -1455, -2307, -980, -2222, 318, 4695, -3133, -1648, -4522, -7550, 4450, 1006, 2052, 5868, 3052, 5824, 6688, -4258, 1356, 7251, -815, 285, 5834, -7636, 7605, -262, 4932, -1902, 4885, -4895, -7056, -628, -3422, 1162, -2516, -4308, 2649, 1318, 7084, -5507, -6927, 2020, 3082, 469, 2311, 4171, 6731, -6691, -6080, -7146, -1078, -98, 4839, 6374, 7535, -3003, 4232, 2181, -5965, -4373, 5866, 4690, -7612, -4581, -5345, -5466, -644, 3498, 6280, 5562, 3741, 2867, 6245, -2767, 3659, -445, -755, 2620, -3237, -4295, 7211, -3104, 2850, -7599, -811, -5561, -7652, 6908, 6920, 3046, 2579, -6396, -2956, 3135, 2730, 7191, 6850, -7192, 2882, 3503, -1801, -5436, -3204, -5159, 5301, -202, -2764, -867, 608, 3180, 6859, 1119, -863, -1583, 1305] );
    
    macro_rules! test_barrett_mul{
        ($func_name:ident, $q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;
                    let qq = $q as i16;

                    let start = -(qq/2+3);
                    let end = -start;

                    let start2 = -(2*qq);
                    let end2 = -start2;

                    for a_base in ( start ..end ).step_by(16){
                        let a_val = a_base as i16;
                        let ar_overq_val = ((a_val as f64 * 65536.0) / q as f64).round() as i16;
                        
                        let a_vec = _mm256_set1_epi16(a_val);
                        let ar_overq_vec = _mm256_set1_epi16(ar_overq_val);

                        for b_start in ( start2 ..end2 ).step_by(16){
                            let mut b_vals = [0i16; 16];
                            for j in 0..16 {
                                b_vals[j] = (b_start as i16).wrapping_add(j as i16);
                            }

                            let b_vec = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);
                            let res_vec = super::$func_name(b_vec, a_vec, ar_overq_vec);
                            let results = dump_m256i(res_vec);

                            for i in 0..16 {
                                let a = a_val as i32;
                                let b = b_vals[i] as i32;
                                
                                let expected = (a * b).rem_euclid(q);
                                let actual = results[i] as i32;

                                assert!((actual - expected) % q == 0, 
                                    "\nFAILED CONGRUENCE: q={}, a={}, b={}\nGot: {}, Expected (mod q): {}\nBarrett Reduction is mathematically incorrect!", 
                                    q, a, b, actual, expected);
                                
                                assert!(actual > -2 * q && actual < 2 * q,
                                    "\nFAILED RANGE: q={}, a={}, b={}\nGot: {}, which is outside (-2q, 2q)", 
                                    q, a, b, actual);
                            }
                        }
                    }
                    println!("SUCCESS: {} passed exhaustive sampling test.", stringify!($func_name));
                }
            }
        };
    }

    // test_barrett_mul!(barrett_mul_7681, 7681);
    // test_barrett_mul!(barrett_mul_10753, 10753);
    // test_barrett_mul!(barrett_mul_11777, 11777);
    // test_barrett_mul!(barrett_mul_12289, 12289);
    // test_barrett_mul!(barrett_mul_13313, 13313);
    // test_barrett_mul!(barrett_mul_15361, 15361);

    macro_rules! test_mont {
        ($func_name:ident, $q:expr, $r_mod_q:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;
                    let q16 = $q as i16;
                    let r_mod_q = $r_mod_q as i32;

                    let startx = -32768;
                    let endx = 32767;

                    let starty = -32768;
                    let endy = 32767;

                    for x_val in ( startx ..endx ).step_by(16){
                        let x_vec = _mm256_set1_epi16(x_val);

                        for y_val in ( starty ..endy ).step_by(16) {
                            let y_vec = _mm256_set1_epi16(y_val);
                            
                            let res_vec = super::$func_name(x_vec, y_vec);
                            let results = dump_m256i(res_vec);

                            for i in 0..16 {
                                let x = x_val as i32;
                                let y = y_val as i32;
                                let actual = results[i] as i32;

                                let left_hand_side = (actual * r_mod_q).rem_euclid(q);
                                let right_hand_side = (x * y).rem_euclid(q);

                                assert_eq!(
                                    left_hand_side, 
                                    right_hand_side, 
                                    "\nFAILED Identity {}: q={}\nInputs: x={}, y={}\nActual Mont Output: {}\nCheck: ({} * {}) % {} = {} != {} (x*y % q)", 
                                    stringify!($func_name), q, x, y, actual, actual, r_mod_q, q, left_hand_side, right_hand_side
                                );

                                assert!(
                                    actual.abs() < q,
                                    "Range error: actual output {} is not in (-3q, 3q) for x={}, y={}", 
                                    actual, x, y
                                );
                            }
                        }
                    }
                    println!("SUCCESS: Exhaustive identity test for {} passed.", stringify!($func_name));
                }
            }
        };
    }

    test_mont!(montproduct_7681, 7681, 4088);
    test_mont!(montproduct_10753, 10753, 1018);
    test_mont!(montproduct_11777, 11777, 6651);
    test_mont!(montproduct_12289, 12289, 4091);
    test_mont!(montproduct_13313, 13313, 12284);
    test_mont!(montproduct_15361, 15361, 4092);
    
    
    macro_rules! test_ibutterfly {
        ($func_name:ident, $q:expr, $c_list:expr) => {
            #[test]
            fn $func_name() {
                if !is_x86_feature_detected!("avx2") { return; }

                unsafe {
                    let q = $q as i32;
                    let c_test_cases: &[i16] = &$c_list;

                    for &c_i16 in c_test_cases {
                        let c = c_i16 as i32;
                        let cr_val = ((c_i16 as f64 * 65536.0) / q as f64).round() as i16;
                        
                        let c_vec = _mm256_set1_epi16(c_i16);
                        let cr_vec = _mm256_set1_epi16(cr_val);

                        let start = -q;
                        let end = -start;

                        for a_base in (start..end).step_by(256) {
                            for b_base in (start..end).step_by(256) {
                                let mut a_vals = [0i16; 16];
                                let mut b_vals = [0i16; 16];
                                
                                for i in 0..16 {
                                    a_vals[i] = (a_base as i16).wrapping_add(i as i16);
                                    b_vals[i] = (b_base as i16).wrapping_add((i * 7) as i16);
                                }

                                let a_vec = _mm256_loadu_si256(a_vals.as_ptr() as *const __m256i);
                                let b_vec = _mm256_loadu_si256(b_vals.as_ptr() as *const __m256i);

                                let [ans0_vec, ans1_vec] = super::$func_name(a_vec, b_vec, c_vec, cr_vec);
                                
                                let ans0_res = dump_m256i(ans0_vec);
                                let ans1_res = dump_m256i(ans1_vec);

                                for i in 0..16 {
                                    let a = a_vals[i] as i32;
                                    let b = b_vals[i] as i32;

                                    let expected0 = (a + b).rem_euclid(q);
                                    let expected1 = ((a - b) * c).rem_euclid(q);

                                    let res0 = (ans0_res[i] as i32).rem_euclid(q);
                                    let res1 = (ans1_res[i] as i32).rem_euclid(q);

                                    assert_eq!(res0, expected0, 
                                        "\nFAILED iNTT ans0: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a, b, ans0_res[i], expected0);

                                    assert_eq!(res1, expected1, 
                                        "\nFAILED iNTT ans1: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a, b, ans1_res[i], expected1);

                                    assert!(res0 > -q && res0 < q,
                                        "\nFAILED iNTT ans1: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a, b, ans1_res[i], expected1);

                                    assert!(res1 > -q && res1 < q,
                                        "\nFAILED iNTT ans1: q={}, c={}\na={}, b={}\nGot: {}, Expected: {}", 
                                        q, c, a, b, ans1_res[i], expected1);
                                }
                            }
                        }
                    }
                    println!("SUCCESS: {} passed for all c cases and i16 range.", stringify!($func_name));
                }
            }
        };
    }

//     test_ibutterfly!(barrett_ibutterfly_7681, 7681, [3383, 1925, 1213, 1213, 1925, 583, 527, 1728, 849, 527, 583, 849, 1728, 2381, 2784, 2446, 1366, 97, -2138, 2132, 2648, 2784, 2381, 1366, 2446, 2138, -97, 2648, 2132, 2381, -2381, 2446, -2446, 97, -97, 2132, -2132, 2784, -2784, 1366, -1366, 2138, -2138, 2648, -2648, 878, 2273, 330, -2645, -2753, -3654, 365, -1846, 1286, 3092, 675, -2268, -1112, -1794, 2399, -3000, 1286, -1286, 3092, -3092, 675, -675, 2268, -2268, 1794, -1794, 1112, -1112, 2399, -2399, 3000, -3000, 878, -878, 2273, -2273, 330, -330, 2645, -2645, 3654, -3654, 2753, -2753, 365, -365, 1846, -1846, -799, -695, -1381, 1875, 2724, 1908, 1382, -2423, -2469, -3380, 693, -1714, -3080, 3477, -732, 3074, -1740, -2774, -584, 1655, 2941, 2508, -3449, 528, -1080, 2516, -3411, -2551, -202, 243, 2881, -766, 1740, -1740, 2774, -2774, 584, -584, 1655, -1655, 2941, -2941, 2508, -2508, 528, -528, 3449, -3449, 2551, -2551, 3411, -3411, 2516, -2516, 1080, -1080, 202, -202, 243, -243, 766, -766, 2881, -2881, 799, -799, 695, -695, 1381, -1381, 1875, -1875, 2724, -2724, 1908, -1908, 2423, -2423, 1382, -1382, 1714, -1714, 693, -693, 3380, -3380, 2469, -2469, 3080, -3080, 3477, -3477, 3074, -3074, 732, -732, -1587, 198, 2900, -2063, -880, -3188, 219, -3501, 405, -2897, 319, 3837, 1633, 1800, 1996, -869, -1155, 2264, -3566, 3073, -1886, 2573, -1220, -2563, -257, 1478, 3141, 3180, -1402, 3789, 3125, -2819, 1125, -3780, -417, 2593, -2990, -707, 1438, -2681, 550, -1848, 1228, 1097, -1591, 2028, 1952, 2044, -2880, -3532, -1682, 1415, 2562, -3078, 648, -3099, 1003, -1853, 2844, 3041, 993, -2722, 1044, -1408, 3780, -3780, 1125, -1125, 2593, -2593, 417, -417, 707, -707, 2990, -2990, 1438, -1438, 2681, -2681, 550, -550, 1848, -1848, 1228, -1228, 1097, -1097, 1952, -1952, 2044, -2044, 2028, -2028, 1591, -1591, 3099, -3099, 648, -648, 3078, -3078, 2562, -2562, 1682, -1682, 1415, -1415, 3532, -3532, 2880, -2880, 1853, -1853, 1003, -1003, 2844, -2844, 3041, -3041, 1408, -1408, 1044, -1044, 993, -993, 2722, -2722, 198, -198, 1587, -1587, 2063, -2063, 2900, -2900, 3188, -3188, 880, -880, 219, -219, 3501, -3501, 405, -405, 2897, -2897, 319, -319, 3837, -3837, 869, -869, 1996, -1996, 1633, -1633, 1800, -1800, 2563, -2563, 1220, -1220, 1886, -1886, 2573, -2573, 3073, -3073, 3566, -3566, 2264, -2264, 1155, -1155, 1478, -1478, 257, -257, 3141, -3141, 3180, -3180, 3125, -3125, 2819, -2819, 3789, -3789, 1402, -1402, 763, -413, -1704, -3799, -2689, 2583, -2668, -669, 2358, -3445, 2922, 321, 1656, -2799, 3694, 185, -1950, 1129, -398, 2259, -3546, -1604, -2359, 62, 1607, 1667, -1968, 1683, 201, -3626, -1979, 2875, 3081, 94, -1193, -3394, 1131, 1035, 3452, -2996, -2173, 542, -3120, -1266, -1437, 702, 1065, 506, -1230, -2012, -1876, -2002, 2197, 2757, 3006, 346, 3139, -3586, -2169, -2372, 296, 2838, 1959, -1406, -3250, -3239, -1897, 3765, -2457, -1189, -1771, 113, 335, 3483, -738, 329, 218, -118, 3280, 2805, -2840, 1211, 3832, -1872, -639, -3376, -674, 1115, 3751, -621, 535, -2811, 3016, 2760, -2252, 1036, 2110, 2481, 1499, -1657, 2395, -1170, -1775, -1717, 3193, -2433, 1885, 1725, -2717, -2546, 572, 536, -217, -3265, -2951, -2067, 3615, 1393, 856, 111, 2050, 793, -1994, 1784, -1459, -3086, -2671, 3137] );
//     test_ibutterfly!(barrett_ibutterfly_10753, 10753, [4489, 67, 321, 321, 67, 3422, 1656, 4679, 3461, 1656, 3422, 3461, 4679, 1560, 2640, 2637, -1154, 4631, -4832, 3010, 2047, 2640, 1560, 1154, -2637, 4832, -4631, 2047, 3010, 1560, -1560, 2637, -2637, 4631, -4631, 3010, -3010, 2640, -2640, 1154, -1154, 4832, -4832, 2047, -2047, -4611, -746, -3783, -2900, -4351, 4191, 1219, 1186, -597, -2436, -1917, -3013, 644, -1641, 2417, 136, 597, -597, 2436, -2436, 3013, -3013, 1917, -1917, 1641, -1641, 644, -644, 2417, -2417, 136, -136, 4611, -4611, 746, -746, 2900, -2900, 3783, -3783, 4191, -4191, 4351, -4351, 1219, -1219, 1186, -1186, -3907, -380, -3954, 3697, 5147, 3314, 753, 3775, -3982, -3712, -1385, 2031, 2603, -3644, 3171, 2353, -1047, -922, 5122, -2744, 94, 2599, 4455, 2085, -559, 3902, -5194, -3362, 946, -841, 1136, 2582, 922, -922, 1047, -1047, 2744, -2744, 5122, -5122, 2085, -2085, 4455, -4455, 2599, -2599, 94, -94, 3362, -3362, 5194, -5194, 559, -559, 3902, -3902, 946, -946, 841, -841, 1136, -1136, 2582, -2582, 380, -380, 3907, -3907, 3697, -3697, 3954, -3954, 3775, -3775, 753, -753, 3314, -3314, 5147, -5147, 1385, -1385, 2031, -2031, 3712, -3712, 3982, -3982, 2603, -2603, 3644, -3644, 3171, -3171, 2353, -2353, -685, 393, -2883, 4825, -5125, -5295, 721, -84, -2004, 4305, -1896, 5232, -2726, -100, 159, -4053, 267, -4980, 3617, -317, -671, -1279, 331, 1945, 3719, 4818, 1854, 216, 118, -2805, -5134, 2847, 2135, -3092, 3256, -2857, 4683, -128, -2177, -1924, 4193, -4627, -1353, -1828, -5178, 3944, -4577, 2830, -5263, -1266, 1202, 2228, -1289, 1207, -339, 5155, 1943, 1444, 1145, 29, 5012, -3592, -2461, 4098, 2857, -2857, 3256, -3256, 3092, -3092, 2135, -2135, 1924, -1924, 2177, -2177, 4683, -4683, 128, -128, 4577, -4577, 2830, -2830, 5178, -5178, 3944, -3944, 1353, -1353, 1828, -1828, 4193, -4193, 4627, -4627, 5155, -5155, 339, -339, 1207, -1207, 1289, -1289, 1266, -1266, 5263, -5263, 1202, -1202, 2228, -2228, 1943, -1943, 1444, -1444, 29, -29, 1145, -1145, 5012, -5012, 3592, -3592, 2461, -2461, 4098, -4098, 4825, -4825, 2883, -2883, 393, -393, 685, -685, 721, -721, 84, -84, 5295, -5295, 5125, -5125, 159, -159, 4053, -4053, 2726, -2726, 100, -100, 1896, -1896, 5232, -5232, 2004, -2004, 4305, -4305, 1279, -1279, 671, -671, 331, -331, 1945, -1945, 317, -317, 3617, -3617, 4980, -4980, 267, -267, 3719, -3719, 4818, -4818, 216, -216, 1854, -1854, 118, -118, 2805, -2805, 5134, -5134, 2847, -2847, -1102, -498, 1437, -1107, 3259, -5182, -3293, -3098, -2129, 2336, 4783, 2854, 5096, -4313, -2664, -1360, -1036, -5308, 4894, -787, 4847, -4864, -2159, 3298, 670, 3210, 1878, 10, -2351, 4946, -1961, 3778, 1196, 3097, -3192, -4861, -4524, -4181, -2024, 549, -4711, 3472, -3800, -3942, 3223, -5262, -881, 2295, -2545, 4819, 283, 1533, -940, -4484, -656, 1538, -3992, -5163, -1825, 1361, 2343, 1293, 607, -4314, -1907, -1135, -774, 1267, 1317, -2137, -3390, -2215, 3661, -3645, -2032, -3104, 2076, 3687, 290, -697, -2758, 3959, 3572, -1985, 1082, 3258, 3226, 2777, -2266, 264, -1280, 3818, -3689, -301, 156, -1339, -4209, -1180, 2425, -3789, 2966, -2160, -4931, 5168, 2037, -4043, -2056, -3310, -2670, -3965, 3911, 3170, -1445, 2546, 38, 1466, -1590, 2482, -4999, -1000, -3429, -5238, 3903, -3930, 815, 2515, 3543, -840] );
//     test_ibutterfly!(barrett_ibutterfly_11777, 11777, [5322, 4976, 4201, 4201, 4976, 2497, -3410, 4578, -337, 3410, -2497, 337, -4578, 1224, -4782, 1447, -293, 1915, 2380, 4525, 5692, 293, -1447, 4782, -1224, 5692, 4525, 2380, 1915, 1224, -1224, 1447, -1447, 1915, -1915, 4525, -4525, 293, -293, 4782, -4782, 5692, -5692, 2380, -2380, -4633, 4212, -4148, 5558, 5573, -5020, 503, -3587, 1265, -4114, 2838, -5722, 3114, -2469, 2333, -3268, 5722, -5722, 2838, -2838, 4114, -4114, 1265, -1265, 2333, -2333, 3268, -3268, 3114, -3114, 2469, -2469, 5558, -5558, 4148, -4148, 4212, -4212, 4633, -4633, 503, -503, 3587, -3587, 5573, -5573, 5020, -5020, 1345, -2326, 3384, -2615, -5739, -5197, -2020, 1961, -1548, -5444, -2244, -690, 3488, -2584, 2500, -3010, -309, 4282, -5203, -2639, 5537, -1860, -1362, 5709, 3238, -2885, -403, 1352, -4055, -5246, -5513, -3679, 5709, -5709, 1362, -1362, 5537, -5537, 1860, -1860, 2639, -2639, 5203, -5203, 309, -309, 4282, -4282, 4055, -4055, 5246, -5246, 3679, -3679, 5513, -5513, 3238, -3238, 2885, -2885, 1352, -1352, 403, -403, 2020, -2020, 1961, -1961, 5197, -5197, 5739, -5739, 3384, -3384, 2615, -2615, 2326, -2326, 1345, -1345, 3488, -3488, 2584, -2584, 3010, -3010, 2500, -2500, 1548, -1548, 5444, -5444, 690, -690, 2244, -2244, -1386, 3890, 4586, -4748, 1596, -2695, -3683, -3998, -1120, -1478, 2599, -5680, -3452, -576, 5491, -4365, 2873, 3560, -1952, 1230, 1534, 2487, 1688, 2315, 4697, -5137, 5622, 5073, -50, 4771, -1483, -1936, -1329, -5039, 5570, -831, -4547, 2601, -2255, -347, 4302, -756, -4993, -3834, -1203, 4322, -1470, 3412, 4779, -4482, 3206, 2541, -2926, -2978, -3062, 3404, 2802, 2562, 1216, -5798, 3673, 2114, 1056, 2403, 1470, -1470, 3412, -3412, 4322, -4322, 1203, -1203, 4302, -4302, 756, -756, 3834, -3834, 4993, -4993, 347, -347, 2255, -2255, 2601, -2601, 4547, -4547, 5039, -5039, 1329, -1329, 5570, -5570, 831, -831, 2562, -2562, 2802, -2802, 5798, -5798, 1216, -1216, 2403, -2403, 1056, -1056, 2114, -2114, 3673, -3673, 4779, -4779, 4482, -4482, 2541, -2541, 3206, -3206, 3062, -3062, 3404, -3404, 2978, -2978, 2926, -2926, 576, -576, 3452, -3452, 5491, -5491, 4365, -4365, 5680, -5680, 2599, -2599, 1478, -1478, 1120, -1120, 1596, -1596, 2695, -2695, 3998, -3998, 3683, -3683, 4748, -4748, 4586, -4586, 1386, -1386, 3890, -3890, 4697, -4697, 5137, -5137, 5073, -5073, 5622, -5622, 1936, -1936, 1483, -1483, 50, -50, 4771, -4771, 3560, -3560, 2873, -2873, 1230, -1230, 1952, -1952, 2315, -2315, 1688, -1688, 2487, -2487, 1534, -1534, 4905, 5199, -3845, 5336, 295, -3649, -2710, 4205, 2550, 3996, 4520, -4971, -371, -4074, -2893, 4007, -1148, 2607, 603, 5822, 3002, 4745, 1765, 4716, -5169, 1654, 1819, 24, 599, 3689, -3879, 1043, -2105, -2883, -4727, -1422, 5857, -2765, 3104, -3643, 3225, 4361, 4451, -4675, 2643, -4308, -2468, 3341, -3131, 1273, 1115, 1578, 5029, 4783, -1109, 1821, 3875, 1223, -3051, 3061, -4819, -3588, -1372, 44, 3435, 3166, 3610, -4133, 3539, 3135, 4765, -3449, 2265, 5318, 51, 551, 2745, 5410, 2200, -2062, 5713, 3628, -1213, 1790, 2172, -5630, -3414, -2603, 5495, -2099, 3074, -1575, 743, 2826, -438, 810, -5518, 5042, 5381, 3982, -1200, 3266, 644, 261, 803, -1485, -3325, -5181, 260, 5811, -1710, 3001, -4071, -3818, -856, -2067, 5824, -1736, -2973, 5795, -3326, 141, 5004, -3491, 1233, -2237, 409, 2047] );
//     test_ibutterfly!(barrett_ibutterfly_12289, 12289, [1479, 4043, 5146, 5146, 4043, 5736, -722, 4134, -1305, 722, -5736, 1305, -4134, 1646, 3621, 1212, 2545, 5860, -3542, 3195, 3504, 2545, 1212, 3621, 1646, 3504, 3195, 3542, -5860, 1646, -1646, 1212, -1212, 5860, -5860, 3195, -3195, 2545, -2545, 3621, -3621, 3504, -3504, 3542, -3542, -3006, 2744, 2975, 563, 2639, 4821, -949, -2625, 5023, 5828, -4591, 5728, -3328, -5777, 1351, 4978, 4591, -4591, 5728, -5728, 5023, -5023, 5828, -5828, 4978, -4978, 1351, -1351, 3328, -3328, 5777, -5777, 2975, -2975, 563, -563, 3006, -3006, 2744, -2744, 949, -949, 2625, -2625, 4821, -4821, 2639, -2639, 4805, -3553, 2294, -1062, 3712, -3135, -2747, -4846, -955, -790, -1170, 2319, 5086, -1326, -3201, 3014, 2963, -4896, -3051, -2366, -4320, 1000, 3091, -81, -1177, 4255, 1635, -2768, -4611, 726, -1853, -140, 1000, -1000, 4320, -4320, 3091, -3091, 81, -81, 2963, -2963, 4896, -4896, 3051, -3051, 2366, -2366, 1853, -1853, 140, -140, 4611, -4611, 726, -726, 4255, -4255, 1177, -1177, 2768, -2768, 1635, -1635, 3712, -3712, 3135, -3135, 2747, -2747, 4846, -4846, 3553, -3553, 4805, -4805, 2294, -2294, 1062, -1062, 1326, -1326, 5086, -5086, 3014, -3014, 3201, -3201, 1170, -1170, 2319, -2319, 955, -955, 790, -790, -3637, -3459, -145, -5542, -5911, 4890, 2731, 3932, 2548, -4231, 3382, -355, -1759, -3707, 5179, -3694, -2426, -334, 1428, -1696, -4632, -5755, 1260, -4388, -2881, -3284, -5092, -2089, 2013, 3289, -729, -3241, 4278, -1673, -5331, 4989, -4177, -3584, -2525, 1381, -2468, -339, -544, -5791, -2842, 480, 1022, 9, -953, -3748, -5767, -827, -118, -2476, 2197, 5067, -4354, -130, 5374, -2837, -4452, -2396, -3949, 3296, 544, -544, 5791, -5791, 339, -339, 2468, -2468, 2842, -2842, 480, -480, 9, -9, 1022, -1022, 4278, -4278, 1673, -1673, 4989, -4989, 5331, -5331, 3584, -3584, 4177, -4177, 1381, -1381, 2525, -2525, 2396, -2396, 4452, -4452, 3296, -3296, 3949, -3949, 130, -130, 4354, -4354, 5374, -5374, 2837, -2837, 5767, -5767, 827, -827, 3748, -3748, 953, -953, 5067, -5067, 2197, -2197, 118, -118, 2476, -2476, 2548, -2548, 4231, -4231, 355, -355, 3382, -3382, 3707, -3707, 1759, -1759, 3694, -3694, 5179, -5179, 5542, -5542, 145, -145, 3637, -3637, 3459, -3459, 5911, -5911, 4890, -4890, 3932, -3932, 2731, -2731, 2089, -2089, 5092, -5092, 2881, -2881, 3284, -3284, 729, -729, 3241, -3241, 3289, -3289, 2013, -2013, 5755, -5755, 4632, -4632, 1260, -1260, 4388, -4388, 334, -334, 2426, -2426, 1696, -1696, 1428, -1428, -2143, -1065, -4645, -404, -5277, -1168, 1207, -3248, -1378, -1912, -435, -4337, -4096, -493, 2381, 5444, 2187, -2566, 6039, 2422, 2437, -3646, 6022, 2987, 4976, 1607, -3780, -875, 4284, 5088, -1002, 5011, -4861, -354, 5698, 2912, 5012, -2481, -1045, 2859, 773, -390, 3833, 3778, -2401, -442, 1067, 5101, -242, -1537, -4143, -4714, -2678, 3704, 5019, -545, 1017, -4885, -5084, 1632, 3066, -27, 1440, -3763, 4057, -3271, 3364, 1689, 2847, 4414, 4372, -2174, -5042, -2305, 4053, 2645, 4895, -1484, 5195, 2780, -3636, 4938, -2704, -5291, -1777, -1663, -4654, -1426, -2166, -3915, -113, -4919, 160, 3149, -3, -4437, -3247, 2686, 3978, -2969, -3510, -5332, -2865, -2370, -3186, 5407, -2126, 1630, -1153, 2884, 4048, -2249, 5559, 420, -1544, 2178, -476, -3531, 4905, 3985, 671, 3000, -243, -3016, 3400, -2399, 3136, -5191] );
//     test_ibutterfly!(barrett_ibutterfly_13313, 13313, [258, 1455, -2626, 2626, -1455, 5275, -6630, 3024, 6476, 6476, 3024, 6630, -5275, 1415, 4215, 5619, 4196, 4690, 5487, 1463, 4468, 5487, 4690, 4468, 1463, 4196, 5619, 4215, 1415, 1415, -1415, 5619, -5619, 4690, -4690, 1463, -1463, 5487, -5487, 4468, -4468, 4196, -4196, 4215, -4215, -3101, 1278, -4330, -1152, 6072, -4358, 5072, 3902, 1173, -3565, 4995, 2651, -5894, -2970, -5375, 2198, 5375, -5375, 2198, -2198, 5894, -5894, 2970, -2970, 1173, -1173, 3565, -3565, 4995, -4995, 2651, -2651, 6072, -6072, 4358, -4358, 5072, -5072, 3902, -3902, 4330, -4330, 1152, -1152, 3101, -3101, 1278, -1278, -772, 519, 3696, -4968, -6168, 6216, 4753, 1478, 1229, -2430, -4253, 5608, -2169, -456, 714, -2170, 3813, -1408, 3606, -1562, -2627, -1197, 2368, -1454, 272, -3611, 3630, 4630, 6105, -4156, 2878, -3004, 4156, -4156, 6105, -6105, 2878, -2878, 3004, -3004, 3611, -3611, 272, -272, 4630, -4630, 3630, -3630, 3813, -3813, 1408, -1408, 1562, -1562, 3606, -3606, 2627, -2627, 1197, -1197, 2368, -2368, 1454, -1454, 1229, -1229, 2430, -2430, 5608, -5608, 4253, -4253, 2169, -2169, 456, -456, 714, -714, 2170, -2170, 6216, -6216, 6168, -6168, 1478, -1478, 4753, -4753, 519, -519, 772, -772, 3696, -3696, 4968, -4968, 3336, 4657, 5365, -382, -4787, -3063, 2386, 3190, 5893, 2712, 743, 5312, 5303, 3065, 5675, -280, 3042, 631, 6194, -492, 275, -4385, -3248, -735, 3904, -4556, 4331, -894, 2935, 1611, -3048, -917, 3523, 3650, -460, 1137, -3152, -1123, -6488, 3534, -225, 4798, -5078, 5450, 694, 5983, -2018, -1437, 2401, 6253, 5336, 5449, 5039, -4612, 3718, 708, 3345, -2335, 2610, 5583, -1969, -2108, 2600, 5150, 5583, -5583, 2610, -2610, 3345, -3345, 2335, -2335, 1969, -1969, 2108, -2108, 5150, -5150, 2600, -2600, 5449, -5449, 5336, -5336, 2401, -2401, 6253, -6253, 3718, -3718, 708, -708, 5039, -5039, 4612, -4612, 3523, -3523, 3650, -3650, 1137, -1137, 460, -460, 3534, -3534, 6488, -6488, 3152, -3152, 1123, -1123, 4798, -4798, 225, -225, 5450, -5450, 5078, -5078, 694, -694, 5983, -5983, 1437, -1437, 2018, -2018, 3042, -3042, 631, -631, 492, -492, 6194, -6194, 3248, -3248, 735, -735, 4385, -4385, 275, -275, 4556, -4556, 3904, -3904, 894, -894, 4331, -4331, 2935, -2935, 1611, -1611, 917, -917, 3048, -3048, 743, -743, 5312, -5312, 2712, -2712, 5893, -5893, 5675, -5675, 280, -280, 5303, -5303, 3065, -3065, 5365, -5365, 382, -382, 4657, -4657, 3336, -3336, 4787, -4787, 3063, -3063, 3190, -3190, 2386, -2386, -2518, -2693, -4293, 2615, -604, 3924, 162, 1857, 2909, 4994, 2628, 939, 3077, 4914, 3867, -789, -407, -1498, 6413, 3742, -3532, -5972, -242, -4129, 5405, 3375, -1872, -3708, 3644, -5071, -3446, 2903, 5730, 597, 3290, 3212, -5240, 6006, -5402, -4149, -4200, 5247, -6036, -333, -4796, -741, -198, 2168, -1600, 97, -1775, 5308, -5781, 442, -2461, 4086, 2288, -4532, -4125, -790, 5691, -3848, -281, -5933, 4282, 223, 174, 4953, -4781, -4611, -6344, 747, -5857, -6576, 1615, -3967, -3798, -5278, -1195, 2111, 3870, -15, -4801, 549, 6259, -3949, -5421, 753, 3625, -3340, -2427, -455, -4407, -5401, -4681, 3785, -49, -671, -4456, -4730, 5358, 2188, -5528, -1733, 3489, -5122, 4242, -2770, 6473, 5909, -5924, -2603, -3196, 838, 3943, 5506, -4816, 4419, -534, -4642, 5121, 3231, 4225, 1616, 2885, -1198, 4080, 913] );
//     test_ibutterfly!(barrett_ibutterfly_15361, 15361, [3968, 1319, 4309, 4309, 1319, 179, 3261, 3666, -5686, 3261, 179, 5686, -3666, 3874, 6372, 4329, -110, 5407, -6841, 4341, -2201, 6372, 3874, 110, -4329, 2201, -4341, 6841, -5407, 3874, -3874, 4329, -4329, 5407, -5407, 4341, -4341, 6372, -6372, 110, -110, 2201, -2201, 6841, -6841, 2395, 5099, -5361, 2537, 7237, 6707, -6422, -1403, 186, 720, 442, -2702, -7467, -2313, 2572, 5992, 186, -186, 720, -720, 2702, -2702, 442, -442, 5992, -5992, 2572, -2572, 2313, -2313, 7467, -7467, 2395, -2395, 5099, -5099, 2537, -2537, 5361, -5361, 1403, -1403, 6422, -6422, 6707, -6707, 7237, -7237, -3519, 243, -2064, -2539, -100, 2586, -6349, 792, 2792, 3375, 3992, -3065, -5046, 7145, -7399, 4361, -5989, 885, -3937, -121, -4805, 3239, 1883, -6298, 1733, -5184, 2051, -2962, 6276, -2987, 1535, 7437, 5989, -5989, 885, -885, 3937, -3937, 121, -121, 6298, -6298, 1883, -1883, 4805, -4805, 3239, -3239, 1535, -1535, 7437, -7437, 6276, -6276, 2987, -2987, 2962, -2962, 2051, -2051, 5184, -5184, 1733, -1733, 3519, -3519, 243, -243, 2064, -2064, 2539, -2539, 6349, -6349, 792, -792, 2586, -2586, 100, -100, 4361, -4361, 7399, -7399, 7145, -7145, 5046, -5046, 3992, -3992, 3065, -3065, 2792, -2792, 3375, -3375, -1100, -2276, 6649, -6966, -7374, -2673, 2793, 7343, 6403, -10, -2993, -2171, -1790, -5938, 1888, 4584, -3702, 4420, -1860, 7200, -7592, -2135, -5002, -1524, 2776, 1331, 4435, -5626, 4907, -6772, 7527, 5352, -5757, -1969, 1102, -5149, -1316, -852, 2435, 11, -6713, -1210, 1554, 6511, -1668, -1967, -1536, 3469, 2802, 3028, -6163, 72, 2473, -2815, 5355, 4377, -2430, -4468, -5332, -5279, -1000, 4862, -7441, 2046, 1969, -1969, 5757, -5757, 1102, -1102, 5149, -5149, 852, -852, 1316, -1316, 11, -11, 2435, -2435, 3469, -3469, 1536, -1536, 1668, -1668, 1967, -1967, 1210, -1210, 6713, -6713, 1554, -1554, 6511, -6511, 1000, -1000, 4862, -4862, 7441, -7441, 2046, -2046, 2430, -2430, 4468, -4468, 5279, -5279, 5332, -5332, 4377, -4377, 5355, -5355, 2473, -2473, 2815, -2815, 72, -72, 6163, -6163, 2802, -2802, 3028, -3028, 2276, -2276, 1100, -1100, 6649, -6649, 6966, -6966, 2673, -2673, 7374, -7374, 7343, -7343, 2793, -2793, 5938, -5938, 1790, -1790, 1888, -1888, 4584, -4584, 2993, -2993, 2171, -2171, 10, -10, 6403, -6403, 7527, -7527, 5352, -5352, 6772, -6772, 4907, -4907, 4435, -4435, 5626, -5626, 1331, -1331, 2776, -2776, 7592, -7592, 2135, -2135, 1524, -1524, 5002, -5002, 3702, -3702, 4420, -4420, 1860, -1860, 7200, -7200, 6440, 6784, 273, -7367, -2353, 2784, -685, 817, -1794, 6449, 3763, -692, -2296, -1455, -2307, -980, -2222, 318, 4695, -3133, -1648, -4522, -7550, 4450, 1006, 2052, 5868, 3052, 5824, 6688, -4258, 1356, 7251, -815, 285, 5834, -7636, 7605, -262, 4932, -1902, 4885, -4895, -7056, -628, -3422, 1162, -2516, -4308, 2649, 1318, 7084, -5507, -6927, 2020, 3082, 469, 2311, 4171, 6731, -6691, -6080, -7146, -1078, -98, 4839, 6374, 7535, -3003, 4232, 2181, -5965, -4373, 5866, 4690, -7612, -4581, -5345, -5466, -644, 3498, 6280, 5562, 3741, 2867, 6245, -2767, 3659, -445, -755, 2620, -3237, -4295, 7211, -3104, 2850, -7599, -811, -5561, -7652, 6908, 6920, 3046, 2579, -6396, -2956, 3135, 2730, 7191, 6850, -7192, 2882, 3503, -1801, -5436, -3204, -5159, 5301, -202, -2764, -867, 608, 3180, 6859, 1119, -863, -1583, 1305] );
   
}