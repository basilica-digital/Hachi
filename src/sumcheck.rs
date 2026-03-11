type MBundle = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);
type WBundle = (Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>, Vec<u32>);

const Q: u32 = 4294967197;

#[inline(always)]
pub unsafe fn ext_mul_val(a: &[u32; 4], b: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    extmul(&mut c, a, b);
    c
}

#[inline(always)]
pub unsafe fn eval_at_2(t0: &[u32; 4], t1: &[u32; 4]) -> [u32; 4] {
    let t1_times_2 = extadd(t1, t1);
    ext_sub(&t1_times_2, t0)
}

#[inline(always)]
pub unsafe fn fold_val(t0: &[u32; 4], t1: &[u32; 4], r: &[u32; 4]) -> [u32; 4] {
    let diff = ext_sub(t1, t0);
    let term = ext_mul_val(&diff, r);
    extadd(t0, &term)
}

// c = a - b
#[inline(always)]
pub unsafe fn ext_sub(a: &[u32; 4], b: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    let q = (1u64 << 32) - 99;
    for i in 0..4 {
        if a[i] >= b[i] {
            c[i] = a[i] - b[i];
        } else {
            c[i] = (q - (b[i] as u64 - a[i] as u64)) as u32;
        }
    }
    c
}

// c = a*b
#[inline(always)]
pub unsafe fn extadd(a: &[u32], b: &[u32]) -> [u32; 4]{
    let mut c = [0u32; 4];
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        c[i] = ((a[i] as u64 + b[i] as u64) % q) as u32;
    }
    c
}

// c = a*b
#[inline(always)]
pub unsafe fn extmul(c: &mut [u32], a: &[u32], b: &[u32]){
    // h(x) = x^4 - 2
    let mut tmp = [0u64; 8];
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        for j in 0..4{
            tmp[i+j] = (tmp[i+j] + (a[i] as u64 * b[j] as u64)) % q;
        }
    }
    for i in 0..4{
        c[i] = ((tmp[i] + 2*tmp[i+4]) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn extmla(c: &mut [u32], a: &[u32], b: &[u32]){
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
pub unsafe fn ext_div_2(a: &[u32; 4]) -> [u32; 4] {
    let mut c = [0u32; 4];
    let inv_2 = 2147483599u64;
    let q = 4294967197u64;
    for i in 0..4 {
        c[i] = ((a[i] as u64 * inv_2) % q) as u32;
    }
    c
}

// P(r) = A*r^2 + B*r + C
#[inline(always)]
pub unsafe fn eval_quadratic(p0: &[u32; 4], p1: &[u32; 4], p2: &[u32; 4], r: &[u32; 4]) -> [u32; 4] {
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
pub unsafe fn ext_base_mla(c: &mut [u32], a: &[u32], b: u32){
    let q = (1u64 << 32) - 99;
    for i in 0..4{
        c[i] = (c[i] as u64 + (a[i] as u64 * b as u64)) % q;
    }
}

pub unsafe fn build_eq_tbl(tau: &[[u32; 4]]) -> Vec<[u32; 4]> {
    let mut eq_table = vec![[1u32, 0, 0, 0]];
    let ext_one = [1u32, 0, 0, 0];

    for t_j in tau {
        let one_minus_t_j = ext_sub(&ext_one, t_j);
        let mut next_table = Vec::with_capacity(eq_table.len() * 2);
        for val in &eq_table {
            let mut left = [0u32; 4];
            extmul(&mut left, val, &one_minus_t_j);
            next_table.push(left);
            let mut right = [0u32; 4];
            extmul(&mut right, val, t_j);
            next_table.push(right);
        }
        eq_table = next_table;
    }
    eq_table
}

#[inline(always)]
pub unsafe fn ringmla_unmod(c: &mut [u32], a: &[u32], b: &[u32]){
    let mut tmp = [0u64; 2048];
    let q = (1u64 << 32) - 99;
    for i in 0..1024 {
        for j in 0..1024{
            tmp[i+j] = (tmp[i+j] + (a[i] as u64 * b[j] as u64)) % q; 
        }
    }
    for i in 0..2048{
        c[i] = ((c[i] as u64 + tmp[i]) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn ringmls_unmod(c: &mut [u32], a: &[u32], b: &[u32]){
    let mut tmp = [0u64; 2048];
    let q = (1u64 << 32) - 99;
    for i in 0..1024 {
        for j in 0..1024{
            tmp[i+j] = (tmp[i+j] + (a[i] as u64 * b[j] as u64)) % q; 
        }
    }
    for i in 0..2048{
        c[i] = ((c[i] as u64 + q - tmp[i]) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn ringmlah_unmod(c: &mut [u32], a: &[u32], b: &[u32]){
    let mut tmp = [0u64; 1024];
    let q = (1u64 << 32) - 99;
    for i in 0..1024 {
        for j in 0..1024{
            if (i+j) >= 1024{
                tmp[i+j-1024] = (tmp[i+j-1024] + (a[i] as u64 * b[j] as u64)) % q; 
            }
        }
    }
    for i in 0..1024{
        c[i] = ((c[i] as u64 + tmp[i]) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn ringmlsh_unmod(c: &mut [u32], a: &[u32], b: &[u32]){
    let mut tmp = [0u64; 1024];
    let q = (1u64 << 32) - 99;
    for i in 0..1024 {
        for j in 0..1024{
            if (i+j) >= 1024{
                tmp[i+j-1024] = (tmp[i+j-1024] + (a[i] as u64 * b[j] as u64)) % q;
            }
        }
    }
    for i in 0..1024{
        c[i] = ((c[i] as u64 + q - tmp[i]) % q) as u32;
    }
}

#[inline(always)]
pub unsafe fn ringdmlah_unmod(c: &mut [u32], a: &[u32], b: &[u32]){
    let mut tmp = [0u64; 1024];
    let q = (1u64 << 32) - 99;
    for i in 0..1024 {
        for j in 0..1024{
            if (i+j) >= 1024{
                tmp[i+j-1024] = (tmp[i+j-1024] + (a[i] as u64 * b[j] as u64)) % q;
            }
        }
    }
    for i in 0..1024{
        c[i] = ((c[i] as u64 + 2*tmp[i]) % q) as u32;
    }
}



pub unsafe fn mle_M(i: u32, u: u32, m: MBundle) -> [u32; 4]{
    let zero = [0u32, 0u32, 0u32, 0u32];
    // D
    if i == 0{
        if u < (1<<13){
            m.0[4*i..4*(i+1)]
        }else {
            zero
        }
    }else if i == 1{
        if u < (2<<13) && u >= (1<<13){
            m.0[4*(u-1<<13)..4*(u-1<<13+1)]
        }else{
            zero
        }
    }else if i == 2{
        if u < (3<<13) && u >= (2<<13){
            m.0[4*(u-2<<13)..4*(u-2<<13+1)]
        }else{
            zero
        }
    }else if i == 3{
        if u < (4<<13) && u >= (3<<13){
            m.0[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    // B
    else if i == 4{
        if u < (5<<13) && u >= (4<<13){
            m.1[4*(u-4<<13)..4*(u-4<<13+1)]
        }else{
            zero
        }
    }
    // bG
    else if i == 5{
        if u < 1<<13{
            m.7[4*i..4*(i+1)]
        }else if u < (2<<13) && u >= (1<<13){
            m.10[4*(u-1<<13)..4*(u-1<<13+1)]
        }else if u < (3<<13) && u >= (2<<13){
            m.9[4*(u-2<<13)..4*(u-2<<13+1)]
        }else if u < (4<<13) && u >= (3<<13){
            m.8[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    else if i == 6{
        if u < 1<<13{
            m.8[4*i..4*(i+1)]
        }else if u < (2<<13) && u >= (1<<13){
            m.7[4*(u-1<<13)..4*(u-1<<13+1)]
        }else if u < (3<<13) && u >= (2<<13){
            m.10[4*(u-2<<13)..4*(u-2<<13+1)]
        }else if u < (4<<13) && u >= (3<<13){
            m.9[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    else if i == 7{
        if u < 1<<13{
            m.9[4*i..4*(i+1)]
        }else if u < (2<<13) && u >= (1<<13){
            m.8[4*(u-1<<13)..4*(u-1<<13+1)]
        }else if u < (3<<13) && u >= (2<<13){
            m.7[4*(u-2<<13)..4*(u-2<<13+1)]
        }else if u < (4<<13) && u >= (3<<13){
            m.10[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    else if i == 8{
        if u < 1<<13{
            m.10[4*i..4*(i+1)]
        }else if u < (2<<13) && u >= (1<<13){
            m.9[4*(u-1<<13)..4*(u-1<<13+1)]
        }else if u < (3<<13) && u >= (2<<13){
            m.8[4*(u-2<<13)..4*(u-2<<13+1)]
        }else if u < (4<<13) && u >= (3<<13){
            m.7[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    // cG, aG
    else if i == 9{
        if u >= 5<<13 {
            m.4[4*(u-5<<13)..4*(u-5<<13+1)]
        }else if u < (1<<13){
            m.11[4*u..4*(u+1)]
        }else{
            zero
        }
    }
    else if i == 10{
        if u >= 5<<13 {
            m.5[4*(u-5<<13)..4*(u-5<<13+1)]
        }else if u < (2<<13) && u >= (1<<13){
            m.11[4*(u-1<<13)..4*(u-1<<13+1)]
        }else{
            zero
        }
    }
    else if i == 11{
        if u >= 5<<13 {
            m.6[4*(u-5<<13)..4*(u-5<<13+1)]
        }else if u < (3<<13) && u >= (2<<13){
            m.11[4*(u-2<<13)..4*(u-2<<13+1)]
        }else{
            zero
        }
    }
    else if i == 12{
        if u >= 5<<13 {
            m.7[4*(u-5<<13)..4*(u-5<<13+1)]
        }else if u < (4<<13) && u >= (3<<13){
            m.11[4*(u-3<<13)..4*(u-3<<13+1)]
        }else{
            zero
        }
    }
    // cG, A
    else if i == 13{
        if u >= 5<<13 {
            m.2[4*(u-5<<13)..4*(u-5<<13+1)]
        }else if u < (5<<13) && u >= (4<<13){
            m.11[4*(u-4<<13)..4*(u-4<<13+1)]
        }else{
            zero
        }
    }
}

pub unsafe fn mle_m(i: u32, u: u32, m: MBundle) -> [u32; 4]{
    // inside M
    if (i < 14) && (u < (6<<13)){
        mle_M(i, u, m)
    }
    // -(alpha^1024 + 1)
    else if (i+6<<13) == u { 
        m.12
    }
    // 0
    else{ 
        [0u32, 0u32, 0u32, 0u32]
    }
}

pub unsafe fn mle_W(u: u32, l: u32, w: WBundle) -> u32{
    if u < (1<<13) {
        w.0[u*1024+l]
    }else if (u < (2<<13) && u >= (1<<13)){
        w.1[(u-1<<13)*1024+l]
    }else if (u < (3<<13) && u >= (2<<13)){
        w.2[(u-2<<13)*1024+l]
    }else if (u < (4<<13) && u >= (3<<13)){
        w.3[(u-3<<13)*1024+l]
    }else if (u < (5<<13) && u >= (4<<13)){
        w.4[(u-4<<13)*1024+l]
    }else if (u < (6<<13) && u >= (5<<13)){
        w.5[(u-5<<13)*1024+l]
    }
}

pub unsafe fn mle_w(u: u32, l: u32, w: WBundle) -> u32{
    // inside W
    if u < (6<<13) {
        mle_W(u, l, w)
    }
    // r
    else{ 
        r[u-(6<<13)*1024+l]
    }
}

pub unsafe fn mle_alpha(y: u32, alpha_vec: Vec<u32>) -> u32{
    alpha_vec[y]
}

pub unsafe fn prove_sumcheck(mut W_table: Vec<[u32; 4]>, mut M_table: Vec<[u32; 4]>, mut alpha_vec: Vec<[u32; 4]>, mut target_sum: [u32; 4]) {
    let mut x_len = 1 << 16;
    let mut y_len = 1 << 10;

    // fold M, W, x16
    for round in 0..16 {
        let half_x = x_len / 2;
        let mut p0 = [0u32; 4];
        let mut p1 = [0u32; 4];
        let mut p2 = [0u32; 4];
        
        for i in 0..half_x {
            let m0 = M_table[i];
            let m1 = M_table[half_x + i];
            let m2 = eval_at_2(&m0, &m1);
            
            for j in 0..y_len {
                let w0 = W_table[i * y_len + j];
                let w1 = W_table[(half_x + i) * y_len + j];
                let w2 = eval_at_2(&w0, &w1);
                
                let alpha = alpha_vec[j];
                
                // P(0): w0 * alpha * m0
                let t0 = ext_mul_val(&ext_mul_val(&w0, &alpha), &m0);
                p0 = extadd(&p0, &t0);
                
                // P(1): w1 * alpha * m1
                let t1 = ext_mul_val(&ext_mul_val(&w1, &alpha), &m1);
                p1 = extadd(&p1, &t1);
                
                // P(2): w2 * alpha * m2
                let t2 = ext_mul_val(&ext_mul_val(&w2, &alpha), &m2);
                p2 = extadd(&p2, &t2);
            }
        }
        
        // Verifier check p0 + p1 == target_sum
        let current_sum = extadd(&p0, &p1);
        if current_sum != target_sum {
            panic!("Sumcheck failed");
        }
        
        // get random challenge r
        let r_raw = generate_random_data_30bit(1, 4);
        let mut r = [0u32; 4];
        r.copy_from_slice(&r_raw[0..4]);
        
        // folding M, W
        for i in 0..half_x {
            M_table[i] = fold_val(&M_table[i], &M_table[half_x + i], &r);
            for j in 0..y_len {
                W_table[i * y_len + j] = fold_val(&W_table[i * y_len + j], &W_table[(half_x + i) * y_len + j], &r);
            }
        }
        x_len = half_x;

        target_sum = eval_quadratic(&p0, &p1, &p2, &r);
    }

    // fold y
    let m_scalar = M_table[0]; 
    
    for round in 0..10 {
        let half_y = y_len / 2;
        let mut p0 = [0u32; 4];
        let mut p1 = [0u32; 4];
        let mut p2 = [0u32; 4];
        
        for j in 0..half_y {
            let w0 = W_table[j];
            let w1 = W_table[half_y + j];
            let w2 = eval_at_2(&w0, &w1);
            
            let a0 = alpha_vec[j];
            let a1 = alpha_vec[half_y + j];
            let a2 = eval_at_2(&a0, &a1);
            
            // P(0): w0 * a0 * m_scalar
            let t0 = ext_mul_val(&ext_mul_val(&w0, &a0), &m_scalar);
            p0 = extadd(&p0, &t0);
            
            // P(1): w1 * a1 * m_scalar
            let t1 = ext_mul_val(&ext_mul_val(&w1, &a1), &m_scalar);
            p1 = extadd(&p1, &t1);
            
            // P(2): w2 * a2 * m_scalar
            let t2 = ext_mul_val(&ext_mul_val(&w2, &a2), &m_scalar);
            p2 = extadd(&p2, &t2);
        }

        // Verifier check p0 + p1 == target_sum
        let current_sum = extadd(&p0, &p1);
        if current_sum != target_sum {
            panic!("Sumcheck failed");
        }
        
        // get random challenge r
        let r_raw = generate_random_data_30bit(1, 4);
        let mut r = [0u32; 4];
        r.copy_from_slice(&r_raw[0..4]);
        
        // folding W, alpha
        for j in 0..half_y {
            W_table[j] = fold_val(&W_table[j], &W_table[half_y + j], &r);
            alpha_vec[j] = fold_val(&alpha_vec[j], &alpha_vec[half_y + j], &r);
        }
        y_len = half_y;

        target_sum = eval_quadratic(&p0, &p1, &p2, &r);
    }

    let final_w = W_table[0];
    let final_alpha = alpha_vec[0];
    let final_m = M_table[0];

    let final_eval = ext_mul_val(&ext_mul_val(&final_w, &final_alpha), &final_m);

    if final_eval == target_sum {
        println!("Success");
    } else {
        println!("Failed");
    }
}

pub fn sumcheck(){
    let n = 1<<10;

    // M
    let D = &mut generate_random_data_30bit(1<<13, n);
    let B = &mut generate_random_data_30bit(1<<13, n);
    let A = &mut generate_random_data_30bit(1<<13, n);
    let a0G = &mut generate_random_data_30bit(1<<13, n);
    let a1G = &mut generate_random_data_30bit(1<<13, n);
    let a2G = &mut generate_random_data_30bit(1<<13, n);
    let a3G = &mut generate_random_data_30bit(1<<13, n);
    let b0G = &mut generate_random_data_30bit(1<<13, n);
    let b1G = &mut generate_random_data_30bit(1<<13, n);
    let b2G = &mut generate_random_data_30bit(1<<13, n);
    let b3G = &mut generate_random_data_30bit(1<<13, n);
    let cG = &mut generate_random_data_30bit(1<<13, n);

    // z
    let w0 = &mut generate_random_data_4bit(1<<13, n);
    let w1 = &mut generate_random_data_4bit(1<<13, n);
    let w2 = &mut generate_random_data_4bit(1<<13, n);
    let w3 = &mut generate_random_data_4bit(1<<13, n);
    let t = &mut generate_random_data_4bit(1<<13, n);
    let z = &mut generate_random_data_30bit(1<<13, n);
    //                       
    let wbundle: WBundle = (w0, w1, w2, w3, t, z);
    
    // y
    let v0 = &mut generate_random_data_30bit(1, n);
    let v1 = &mut generate_random_data_30bit(1, n);
    let v2 = &mut generate_random_data_30bit(1, n);
    let v3 = &mut generate_random_data_30bit(1, n);
    let u  = &mut generate_random_data_30bit(1, n);
    let u0 = &mut generate_random_data_30bit(1, n);
    let u1 = &mut generate_random_data_30bit(1, n);
    let u2 = &mut generate_random_data_30bit(1, n);
    let u3 = &mut generate_random_data_30bit(1, n);

    // r
    let mut r = vec![0u32; 14*n];
    
    // compute r: Mz = y + r*(x^1024+1)
    for i in 0..(1<<13){
        let o = i * n;
        // v0 = Dw0
        ringmlah_unmod(&mut r[0..n], &D[o..o+n], &w0[o..o+n]);
        // v1 = Dw1
        ringmlah_unmod(&mut r[n..2*n], &D[o..o+n], &w1[o..o+n]);
        // v2 = Dw2
        ringmlah_unmod(&mut r[2*n..3*n], &D[o..o+n], &w2[o..o+n]);
        // v3 = Dw3
        ringmlah_unmod(&mut r[3*n..4*n], &D[o..o+n], &w3[o..o+n]);
        // u = Bt
        ringmlah_unmod(&mut r[4*n..5*n], &B[o..o+n], &t[o..o+n]);
        // u0 = b0w0 - b3w1 - b2w2 - b1w3
        ringmlah_unmod(&mut r[5*n..6*n], &b0G[o..o+n], &w0[o..o+n]);
        ringdmlah_unmod(&mut r[5*n..6*n], &b3G[o..o+n], &w1[o..o+n]);
        ringdmlah_unmod(&mut r[5*n..6*n], &b2G[o..o+n], &w2[o..o+n]);
        ringdmlah_unmod(&mut r[5*n..6*n], &b1G[o..o+n], &w3[o..o+n]);
        // u1 = b1w0 + b0w1 - b3w2 - b2w3
        ringmlah_unmod(&mut r[6*n..7*n], &b1G[o..o+n], &w0[o..o+n]);
        ringmlah_unmod(&mut r[6*n..7*n], &b0G[o..o+n], &w1[o..o+n]);
        ringdmlah_unmod(&mut r[6*n..7*n], &b3G[o..o+n], &w2[o..o+n]);
        ringdmlah_unmod(&mut r[6*n..7*n], &b2G[o..o+n], &w3[o..o+n]);
        // u2 = b2w0 + b1w1 + b0w2 - b3w3
        ringmlah_unmod(&mut r[7*n..8*n], &b2G[o..o+n], &w0[o..o+n]);
        ringmlah_unmod(&mut r[7*n..8*n], &b1G[o..o+n], &w1[o..o+n]);
        ringmlah_unmod(&mut r[7*n..8*n], &b0G[o..o+n], &w2[o..o+n]);
        ringdmlah_unmod(&mut r[7*n..8*n], &b3G[o..o+n], &w3[o..o+n]);
        // u3 = b3w0 + b2w1 + b1w2 + b0w3
        ringmlah_unmod(&mut r[8*n..9*n], &b3G[o..o+n], &w0[o..o+n]);
        ringmlah_unmod(&mut r[8*n..9*n], &b2G[o..o+n], &w1[o..o+n]);
        ringmlah_unmod(&mut r[8*n..9*n], &b1G[o..o+n], &w2[o..o+n]);
        ringmlah_unmod(&mut r[8*n..9*n], &b0G[o..o+n], &w3[o..o+n]);
        // 0 = cw0 - a0z
        ringmlah_unmod(&mut r[9*n..10*n], &cG[o..o+n], &w0[o..o+n]);
        ringmlsh_unmod(&mut r[9*n..10*n], &a0G[o..o+n], &z[o..o+n]);
        // 0 = cw1 - a1z
        ringmlah_unmod(&mut r[10*n..11*n], &cG[o..o+n], &w1[o..o+n]);
        ringmlsh_unmod(&mut r[10*n..11*n], &a1G[o..o+n], &z[o..o+n]);
        // 0 = cw2 - a2z
        ringmlah_unmod(&mut r[11*n..12*n], &cG[o..o+n], &w2[o..o+n]);
        ringmlsh_unmod(&mut r[11*n..12*n], &a2G[o..o+n], &z[o..o+n]);
        // 0 = cw3 - a3z
        ringmlah_unmod(&mut r[12*n..13*n], &cG[o..o+n], &w3[o..o+n]);
        ringmlsh_unmod(&mut r[12*n..13*n], &a3G[o..o+n], &z[o..o+n]);
        // 0 = ct - Az
        ringmlah_unmod(&mut r[13*n..14*n], &cG[o..o+n], &t[o..o+n]);
        ringmlsh_unmod(&mut r[13*n..14*n], &A[o..o+n], &z[o..o+n]);
    }

    // compute alpha
    let alpha = generate_random_data_30bit(1, 4);
    let mut alpha_vec = vec![0u32; 4*1025];
    alpha_vec[0] = 1;
    alpha_vec[1] = 0;
    alpha_vec[2] = 0;
    alpha_vec[3] = 0;
    for i in 1..1024{
        extmul(&mut alpha_vec[4*i..4*(i+1)], &alpha_vec[4*(i-1)..4*i], &alpha);
    }
    let q = (1u64 << 32) - 99;
    // alpha_vec[4*1025] = - (alpha^1024 + 1)
    alpha_vec[4096] = (q - ((alpha_vec[4092] as u64 + 1u664)%q)) as u32;
    alpha_vec[4097] = (q - (alpha_vec[4093] as u64)) as u32;
    alpha_vec[4098] = (q - (alpha_vec[4094] as u64)) as u32;
    alpha_vec[4099] = (q - (alpha_vec[4095] as u64)) as u32;

    // substitute alpha into M
    // M(alpha)
    let mut D_alpha = vec![0u32; 1<<15];
    let mut B_alpha = vec![0u32; 1<<15];
    let mut A_alpha = vec![0u32; 1<<15];
    let mut a0G_alpha = vec![0u32; 1<<15];
    let mut a1G_alpha = vec![0u32; 1<<15];
    let mut a2G_alpha = vec![0u32; 1<<15];
    let mut a3G_alpha = vec![0u32; 1<<15];
    let mut b0G_alpha = vec![0u32; 1<<15];
    let mut b1G_alpha = vec![0u32; 1<<15];
    let mut b2G_alpha = vec![0u32; 1<<15];
    let mut b3G_alpha = vec![0u32; 1<<15];
    let mut cG_alpha = vec![0u32; 1<<15];
    for i in 0..(1<<13){
        for j in 0..n{
            ext_base_mla(&mut D_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], D[i*n+j]);
            ext_base_mla(&mut B_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], B[i*n+j]);
            ext_base_mla(&mut A_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], A[i*n+j]);
            ext_base_mla(&mut a0G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], a0G[i*n+j]);
            ext_base_mla(&mut a1G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], a1G[i*n+j]);
            ext_base_mla(&mut a2G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], a2G[i*n+j]);
            ext_base_mla(&mut a3G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], a3G[i*n+j]);
            ext_base_mla(&mut b0G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], b0G[i*n+j]);
            ext_base_mla(&mut b1G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], b1G[i*n+j]);
            ext_base_mla(&mut b2G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], b2G[i*n+j]);
            ext_base_mla(&mut b3G_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], b3G[i*n+j]);
            ext_base_mla(&mut cG_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], cG[i*n+j]);
        }
    }
    //                          0       1        2         3          4          5          6          7           8          9         10        11         12        
    let mbundle: MBundle = (D_alpha, B_alpha, A_alpha, a0G_alpha, a1G_alpha, a2G_alpha, a3G_alpha, b0G_alpha, b1G_alpha, b2G_alpha, b3G_alpha, cG_alpha, alpha_vec[4096..4100].to_vec());

    // // z(alpha)
    // let mut w0_alpha = vec![0u32; 1<<15];
    // let mut w1_alpha = vec![0u32; 1<<15];
    // let mut w2_alpha = vec![0u32; 1<<15];
    // let mut w3_alpha = vec![0u32; 1<<15];
    // let mut t_alpha = vec![0u32; 1<<15];
    // let mut z_alpha = vec![0u32; 1<<15];
    // for i in 0..(1<<13){
    //     for j in 0..n{
    //         ext_base_mla(&mut w0_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], w0[i*n+j]);
    //         ext_base_mla(&mut w1_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], w1[i*n+j]);
    //         ext_base_mla(&mut w2_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], w2[i*n+j]);
    //         ext_base_mla(&mut w3_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], w3[i*n+j]);
    //         ext_base_mla(&mut t_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], t[i*n+j]);
    //         ext_base_mla(&mut z_alpha[4*i..4*(i+1)], &alpha_vec[4*j..4*(j+1)], z[i*n+j]);
    //     }
    // }
    
    // y(alpha)
    let mut v0_alpha = vec![0u32; 1<<2];
    let mut v1_alpha = vec![0u32; 1<<2];
    let mut v2_alpha = vec![0u32; 1<<2];
    let mut v3_alpha = vec![0u32; 1<<2];
    let mut u_alpha  = vec![0u32; 1<<2];
    let mut u0_alpha = vec![0u32; 1<<2];
    let mut u1_alpha = vec![0u32; 1<<2];
    let mut u2_alpha = vec![0u32; 1<<2];
    let mut u3_alpha = vec![0u32; 1<<2];

    for j in 0..n{
        ext_base_mla(&mut v0_alpha, &alpha_vec[4*j..4*(j+1)], v0[j]);
        ext_base_mla(&mut v1_alpha, &alpha_vec[4*j..4*(j+1)], v1[j]);
        ext_base_mla(&mut v2_alpha, &alpha_vec[4*j..4*(j+1)], v2[j]);
        ext_base_mla(&mut v3_alpha, &alpha_vec[4*j..4*(j+1)], v3[j]);
        ext_base_mla(&mut u_alpha,  &alpha_vec[4*j..4*(j+1)],  u[j]);
        ext_base_mla(&mut u0_alpha, &alpha_vec[4*j..4*(j+1)], u0[j]);
        ext_base_mla(&mut u1_alpha, &alpha_vec[4*j..4*(j+1)], u1[j]);
        ext_base_mla(&mut u2_alpha, &alpha_vec[4*j..4*(j+1)], u2[j]);
        ext_base_mla(&mut u3_alpha, &alpha_vec[4*j..4*(j+1)], u3[j]);
    }

    // // r(alpha)
    // let mut r_alpha = vec![0u32; 14*4];

    // for i in 0..14{
    //     for j in 0..n{
    //         ext_base_mla(&mut r_alpha[i*4..(i+1)*4], &alpha_vec[4*j..4*(j+1)], r[i*n+j]);
    //     }
    // }

    // MLE
    let tau_0 = [[0u32; 4]; 26]; // log 1024 + log 6*(1<<13)
    let tau_1 = [[0u32; 4]; 4]; // log 14
    for i in 0..26{
        if i < 4{
            tau_1[i] = generate_random_data_30bit(1, 4);
        }
        tau_0[i] = generate_random_data_30bit(1, 4);
    }
    let tbl_tau0 = build_eq_tbl(&tau_0);
    let tbl_tau1 = build_eq_tbl(&tau_1);
    let mut a = [0u32; 4];
    extmla(&mut a, &tbl_tau1[0], &v0_alpha);
    extmla(&mut a, &tbl_tau1[1], &v1_alpha);
    extmla(&mut a, &tbl_tau1[2], &v2_alpha);
    extmla(&mut a, &tbl_tau1[3], &v3_alpha);
    extmla(&mut a, &tbl_tau1[4], &u_alpha);
    extmla(&mut a, &tbl_tau1[5], &u0_alpha);
    extmla(&mut a, &tbl_tau1[6], &u1_alpha);
    extmla(&mut a, &tbl_tau1[7], &u2_alpha);
    extmla(&mut a, &tbl_tau1[8], &u3_alpha);

    let mut M_table = vec![[0u32; 4]; 1<<16];
    for x in 0..(6*(1<<13) + 14) {
        let mut sum = [0u32; 4];
        for i in 0..14 {
            let m_val = mle_m(i, x as u32, wbundle.clone());
            let mut term = [0u32; 4];
            extmul(&mut term, &tbl_tau1[i], &m_val);
            sum = extadd(&sum, &term);
        }
        M_table[x] = sum;
    }

    let mut W_table = vec![[0u32; 4]; 1<<26]; // (6*(1<<13) + 14) * 1024 padding to 1<<26
    for u in 0..(6*(1<<13) + 14){
        for l in 0..1024{
            W_table[u*1024 + l] = [mle_w(u as u32, l as u32, wbundle.clone()), 0, 0, 0];
        }
    }

    // Sumcheck
    prove_sumcheck(W_table, M_table, alpha_vec, a);
}