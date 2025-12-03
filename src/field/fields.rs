#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
// RNS: [257, 3329, 7681, 7937, 9473, 10753]

// barrett_fake
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_257(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(85));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(257));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_3329(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(20));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(3329));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(9));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(7681));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7937(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(8));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(7937));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_9473(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(7));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(9473));
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

""" Sage
def to_int16(x):
    u = Integer(x) & 0xFFFF 
    if u >= 0x8000:
        return u - 0x10000
    return u

f = [769, 3329, 7681, 7937, 9473, 10753]
for i in f:
    for j in range(10000):
        a = ZZ.random_element(-i, i)
        r = ZZ.random_element(-(2^15), 2^15-1)
        m = ZZ.random_element(-(2^15), 2^15-1)
        m1 = round((m*(2^16))/(i))
        t = (r*m1) >> 16
        d = to_int16(r*m)
        dp = to_int16(a+d)
        ds = to_int16(a-d)
        ansp = to_int16(dp - to_int16(t*i))
        anss = to_int16(ds + to_int16(t*i))
        if ansp%i != (a + r*m)%i or ansp - (a + (r*m)%i) < -3*i or ansp - (a + (r*m)%i) > 3*i:
            print(ansp, a + (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
        if anss%i != (a - r*m)%i or anss - (a - (r*m)%i) < -3*i or anss - (a - (r*m)%i) > 3*i:
            print(anss, a - (r*m)%i, "||", a, r, m, ansp - (a + (r*m)%i))
"""

// barrett_butterfly
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_257(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(257));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_3329(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(3329));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(7681));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_7937(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(7937));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_9473(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(9473));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_butterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let t = _mm256_mulhi_epi16(b, cr);
    let d = _mm256_mullo_epi16(b, c);
    let dp = _mm256_add_epi16(a, d);
    let ds = _mm256_sub_epi16(a, d);
    let e = _mm256_sub_epi16(t, _mm256_set1_epi16(10753));
    let ansp = _mm256_sub_epi16(dp, e);
    let anss = _mm256_add_epi16(dp, e);
    [ansp, anss]
}

"""

"""
// barrett_fake_32 
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_257_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(5585133);
    let m2 = _mm256_set1_epi32(257);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(5585133));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(257));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(5585133));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(257));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_3329_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(1290167);
    let m2 = _mm256_set1_epi32(3329);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(1290167));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(3329));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(1290167));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(3329));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7681_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(559168);
    let m2 = _mm256_set1_epi32(7681);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(559168));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(7681));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(559168));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(7681));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_7937_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(541132);
    let m2 = _mm256_set1_epi32(7937);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(541132));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(7937));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(541132));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(7937));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_9473_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(453390);
    let m2 = _mm256_set1_epi32(9473);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(453390));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(9473));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(453390));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(9473));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_10753_32(x: __m256i, y: __m256i) -> __m256i{
    let m1 = _mm256_set1_epi32(399420);
    let m2 = _mm256_set1_epi32(10753);
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(399420));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(10753));
    let f = _mm256_sub_epi16(x, e);
    let d1 = _mm256_mulhi_epi32(y, _mm256_set1_epi32(399420));
    let e1 = _mm256_mullo_epi32(d1, _mm256_set1_epi32(10753));
    let f1 = _mm256_sub_epi16(y, e1);
    let ans = pack_32_16(f, f1);
}

// barrett_mul
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_257(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(257));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_3329(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(3329));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7681(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(7681));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7937(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(7937));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_9473(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(9473));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_10753(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(10753));
    let g = _mm256_sub_epi16(d, f);
    g
}

// barrett_ibutterfly
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_257(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_257(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_3329(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_3329(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_7681(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_7681(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_7937(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_7937(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_9473(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_9473(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_ibutterfly_10753(a: __m256i, b: __m256i, c: __m256i, cr:__m256i) -> [__m256i, __m256i]{
    let ans0 = barrett_fake_10753(_mm256_add_epi16(a, b));
    let tmp = _mm256_sub_epi16(a, b);
    let ans1 = barrett_mul(tmp, c, cr);
    [ans0, ans1]
}
