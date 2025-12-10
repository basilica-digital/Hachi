#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
// RNS: [257, 3329, 7681, 7937, 9473, 10753]

// reduce 2^32-99
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

// montproduct

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_7681(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(57857));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(7681));
    let f = _mm256_sub_epi16(hi, e);
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_10753(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(54785));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(10753));
    let f = _mm256_sub_epi16(hi, e);
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_11777(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(53761));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(11777));
    let f = _mm256_sub_epi16(hi, e);
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_12289(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(53249));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(12289));
    let f = _mm256_sub_epi16(hi, e);
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_13313(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(52225));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(13313));
    let f = _mm256_sub_epi16(hi, e);
}

#[target_feature(enable = "avx2")]
pub unsafe fn montproduct_15361(x: __m256i, y: __m256i) -> __m256i{
    let lo = _mm256_mullo_epi16(x, y);
    let hi = _mm256_mulhi_epi16(x, y);
    let d = _mm256_mullo_epi16(lo, _mm256_set1_epi16(50177));
    let e = _mm256_mulhi_epi16(d, _mm256_set1_epi16(15361));
    let f = _mm256_sub_epi16(hi, e);
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

    #[test]
    fn test_barrett_7681_correctness() {
        if !is_x86_feature_detected!("avx2") {
            println!("Skipping AVX2 test: Instruction set not supported");
            return;
        }

        unsafe {
            let inputs: [i16; 16] = [
                0, 1, 10, 7680,
                7681, 7682, 15362,
                -1, -10, -7681,
                32767, -32768,
                100, 200, 300, 400
            ];

            let input_vec = _mm256_loadu_si256(inputs.as_ptr() as *const __m256i);
            let output_vec = barrett_fake_7681(input_vec);
            let outputs = dump_m256i(output_vec);

            println!("| Index |   Input |  Output | Diff (In - Out) | Is Multiple of 7681? |");
            println!("|-------|---------|---------|-----------------|----------------------|");

            for i in 0..16 {
                let x = inputs[i];
                let y = outputs[i];
                let diff = (x as i32) - (y as i32);
                let is_valid = diff % 7681 == 0;
                println!("| {:5} | {:7} | {:7} | {:15} | {:^20} |", 
                    i, x, y, diff, if is_valid { "Ok" } else { "Failed" });
                assert!(is_valid, "Reduction failed at index {}: input={}, output={}", i, x, y);
            }
        }
    }
}