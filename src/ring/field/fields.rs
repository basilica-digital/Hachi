#[cfg(target_arch = "x86_64")]
use std::arch::x86_64::*;
use crate::field::macros::*;
// RNS fields: [26113, 25601, 23041, 19457, 18433]

// barrett_fake
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_26113(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(3));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(26113));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_25601(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(3));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(25601));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_23041(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(3));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(23041));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_19457(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(4));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(19457));
    let f = _mm256_sub_epi16(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_18433(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi16(x, _mm256_set1_epi16(4));
    let e = _mm256_mullo_epi16(d, _mm256_set1_epi16(18433));
    let f = _mm256_sub_epi16(x, e);
    f
}

// barrett_fake_32
#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_26113_32(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(164476));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(26113));
    let f = _mm256_sub_epi32(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_25601_32(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(167765));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(25601));
    let f = _mm256_sub_epi32(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_23041_32(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(186405));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(23041));
    let f = _mm256_sub_epi32(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_19457_32(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(220741));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(19457));
    let f = _mm256_sub_epi32(x, e);
    f
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_fake_18433_32(x: __m256i) -> __m256i{
    let d = _mm256_mulhi_epi32(x, _mm256_set1_epi32(233004));
    let e = _mm256_mullo_epi32(d, _mm256_set1_epi32(18433));
    let f = _mm256_sub_epi32(x, e);
    f
}


// a*b % r
// #[target_feature(enable = "avx2")]
// pub unsafe fn dot_product_26113(a: __m256i, b: __m256i) -> __m256i{

// }

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_26113(b: __m256i, a: __m256i, ar_overq: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, ar_overq);
    let d = _mm256_mullo_epi16(b, a);
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_8534_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(21418));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(8534));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7720_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(19375));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(7720));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_619_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(1554));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(619));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_5428_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(13623));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(5428));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_1910_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(4794));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(1910));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_8645_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(21696));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(8645));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7205_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(18082));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(7205));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_5700_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(14305));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(5700));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_4719_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(11843));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(4719));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_3045_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(7642));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(3045));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_3595_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(9022));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(3595));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_7249_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(18193));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(7249));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_1269_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(3185));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(1269));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_2121_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(5323));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(2121));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}

#[target_feature(enable = "avx2")]
pub unsafe fn barrett_mul_4305_26113(b: __m256i) -> __m256i{
    let t = _mm256_mulhi_epi16(b, _mm256_set1_epi16(10804));
    let d = _mm256_mullo_epi16(b, _mm256_set1_epi16(4305));
    let f = _mm256_mullo_epi16(t, _mm256_set1_epi16(26113));
    let g = _mm256_sub_epi16(d, f);
    g
}