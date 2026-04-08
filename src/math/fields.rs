use crate::Q;

#[inline(always)]
pub fn ext_mul_val(a: &[u32; 4], b: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    extmul(&mut c, a, b);
    c
}

#[inline(always)]
pub fn fold_val(t0: &[u32; 4], t1: &[u32; 4], r: &[u32; 4]) -> [u32; 4]{
    let diff = ext_sub(t1, t0);
    let term = ext_mul_val(&diff, r);
    extadd(t0, &term)
}

// c = a - b
#[inline(always)]
pub fn ext_sub(a: &[u32; 4], b: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    let q = 4294967197u32;
    for i in 0..4 {
        c[i] = if a[i] >= b[i] { a[i] - b[i] } else { a[i] + q - b[i] };
    }
    c
}

// c = a*b
#[inline(always)]
pub fn extadd(a: &[u32], b: &[u32]) -> [u32; 4] {
    let mut c = [0u32; 4];
    let q = 4294967197u64;
    for i in 0..4 {
        let sum = a[i] as u64 + b[i] as u64;
        c[i] = if sum >= q { (sum - q) as u32 } else { sum as u32 };
    }
    c
}

#[inline(always)]
pub fn extdbl(a: &[u32; 4]) -> [u32; 4]{
    let mut c = [0u32; 4];
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        c[i] = ((a[i] as u64 * 2u64) % q) as u32;
    }
    c
}

#[inline(always)]
pub fn extneg(a: &[u32; 4]) -> [u32; 4]{
    let mut c = [0u32; 4];
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        c[i] = ((q - a[i] as u64) % q) as u32;
    }
    c
}


#[inline(always)]
pub fn extmul(c: &mut [u32], a: &[u32], b: &[u32]) {
    #[inline(always)]
    fn red(x: u64) -> u64 { (x as u32 as u64) + (x >> 32) * 99 }

    let a0 = a[0] as u64; let a1 = a[1] as u64; let a2 = a[2] as u64; let a3 = a[3] as u64;
    let b0 = b[0] as u64; let b1 = b[1] as u64; let b2 = b[2] as u64; let b3 = b[3] as u64;

    let p00 = red(a0 * b0); let p01 = red(a0 * b1); let p02 = red(a0 * b2); let p03 = red(a0 * b3);
    let p10 = red(a1 * b0); let p11 = red(a1 * b1); let p12 = red(a1 * b2); let p13 = red(a1 * b3);
    let p20 = red(a2 * b0); let p21 = red(a2 * b1); let p22 = red(a2 * b2); let p23 = red(a2 * b3);
    let p30 = red(a3 * b0); let p31 = red(a3 * b1); let p32 = red(a3 * b2); let p33 = red(a3 * b3);

    let t0 = p00;
    let t1 = p01 + p10;
    let t2 = p02 + p11 + p20;
    let t3 = p03 + p12 + p21 + p30;
    let t4 = p13 + p22 + p31;
    let t5 = p23 + p32;
    let t6 = p33;

    #[inline(always)]
    fn final_red(mut x: u64) -> u32 {
        x = (x as u32 as u64) + (x >> 32) * 99;
        let mut res = (x as u32 as u64) + (x >> 32) * 99;
        if res >= 4294967197 { res -= 4294967197; }
        res as u32
    }

    c[0] = final_red(t0 + 2 * t4);
    c[1] = final_red(t1 + 2 * t5);
    c[2] = final_red(t2 + 2 * t6);
    c[3] = final_red(t3);
}

#[inline(always)]
pub fn extmla(c: &mut [u32], a: &[u32], b: &[u32]){
    // h(x) = x^4 - 2
    let mut tmp = [0u64; 8];
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        for j in 0..4{
            tmp[i+j] = (tmp[i+j] + (a[i] as u64 * b[j] as u64)) % q;
        }
    }
    for i in 0..4{
        c[i] = ((c[i] as u64 + tmp[i] + 2*tmp[i+4]) % q) as u32;
    }
}

#[inline(always)]
pub fn ext_div_2(a: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    let inv_2 = 2147483599u64;
    let q = 4294967197u64;
    for i in 0..4 {
        c[i] = ((a[i] as u64 * inv_2) % q) as u32;
    }
    c
}

#[inline(always)]
pub fn ext_base_mla(c: &mut [u32], a: &[u32], b: u32){
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        c[i] = ((c[i] as u64 + (a[i] as u64 * b as u64)) % q) as u32;
    }
}

// P(r) = A*r^2 + B*r + C
#[inline(always)]
pub fn eval_quadratic(p0: &[u32; 4], p1: &[u32; 4], p2: &[u32; 4], r: &[u32; 4]) -> [u32; 4] {
    // A = (P(2) - 2*P(1) + P(0)) / 2
    let p1_times_2 = extadd(p1, p1);
    let mut a_coeff = ext_sub(p2, &p1_times_2);
    a_coeff = extadd(&a_coeff, p0);
    a_coeff = ext_div_2(&a_coeff);

    // B = P(1) - P(0) - A
    let mut b_coeff = ext_sub(p1, p0);
    b_coeff = ext_sub(&b_coeff, &a_coeff);

    // P(r) = (A * r + B) * r + P(0)
    let mut res = ext_mul_val(&a_coeff, r);
    res = extadd(&res, &b_coeff);
    res = ext_mul_val(&res, r);
    res = extadd(&res, p0);
    res
}

#[inline(always)]
pub fn ext_scalar_mul(a: &[u32; 4], b: u32) -> [u32; 4] {
    let mut c = [0u32; 4];
    for i in 0..4 {
        c[i] = ((a[i] as u64 * b as u64) % Q as u64) as u32;
    }
    c
}

pub fn eval_degree_18(p: &[[u32; 4]; 19], r: &[u32; 4], w_i: &[u32; 19]) -> [u32; 4] {
    if r[1] == 0 && r[2] == 0 && r[3] == 0 && r[0] < 19 {
        return p[r[0] as usize];
    }
    
    let mut pref = [[0u32; 4]; 19];
    let mut suff = [[0u32; 4]; 19];
    
    let mut current_pref = [1u32, 0, 0, 0];
    for i in 0..19 {
        pref[i] = current_pref;
        let mut r_minus_i = *r;
        r_minus_i[0] = ((r_minus_i[0] as u64 + Q as u64 - i as u64) % Q as u64) as u32;
        current_pref = ext_mul_val(&current_pref, &r_minus_i);
    }
    
    let mut current_suff = [1u32, 0, 0, 0];
    for i in (0..19).rev() {
        suff[i] = current_suff;
        let mut r_minus_i = *r;
        r_minus_i[0] = ((r_minus_i[0] as u64 + Q as u64 - i as u64) % Q as u64) as u32;
        current_suff = ext_mul_val(&current_suff, &r_minus_i);
    }
    
    let mut res = [0u32; 4];
    for i in 0..19 {
        let mut l_i = ext_mul_val(&pref[i], &suff[i]);
        l_i = ext_scalar_mul(&l_i, w_i[i]);
        let term = ext_mul_val(&l_i, &p[i]);
        res = extadd(&res, &term);
    }
    res
}

pub fn mod_inverse(a: u64, m: u64) -> u64 {
    let mut mn = (m as i64, a as i64);
    let mut xy = (0, 1);
    while mn.1 != 0 {
        xy = (xy.1, xy.0 - (mn.0 / mn.1) * xy.1);
        mn = (mn.1, mn.0 % mn.1);
    }
    if xy.0 < 0 {
        (xy.0 + m as i64) as u64
    } else {
        xy.0 as u64
    }
}