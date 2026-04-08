use tfhe_ntt::prime32::Plan;
use crate::utils::random::generate_random_q_element;
use crate::math::fields::mod_inverse;
use crate::Q;
use crate::utils::ds::AlignedU32Vec;

pub struct SetupParams {
    pub constraints: usize,
    pub height: usize,
    pub n: usize,
    pub q: u32,
    pub q64: u64,
    pub plans: Vec<Plan>,
    pub g_matrix: Vec<u32>,
    pub d0: AlignedU32Vec,
    pub d1: AlignedU32Vec,
    pub d2: AlignedU32Vec,
    pub d3: AlignedU32Vec,
    pub b: AlignedU32Vec,
    pub a: AlignedU32Vec,
    pub w_i: [u32; 19],
}

fn compute_w_i() -> [u32; 19] {
    let mut w = [0u32; 19];
    for i in 0..19usize {
        let mut denom = 1u64;
        for j in 0..19usize {
            if i != j {
                let diff = if i > j { 
                    (i - j) as u64 
                } else { 
                    Q as u64 - (j - i) as u64 
                };
                denom = (denom * diff) % Q as u64;
            }
        }
        w[i] = mod_inverse(denom, Q as u64) as u32;
    }
    w
}

pub fn setup() -> SetupParams{
    let constraints = 11;
    let height = 1 << 13;
    let n = 1 << 10;
    let q = 4294967197u32;
    let q64 = 4294967197u64;
    let moduli = [2079301633, 2079305729];
    let plans: Vec<Plan> = moduli.iter().map(|&m| Plan::try_new(1<<11, m).unwrap()).collect();
    let mut g_matrix = vec![0u32; 8];
    for i in 0..8 {
        g_matrix[i] = 1 << (i * 4);
    }
    let d0 = generate_random_q_element(height, n);
    let d1 = generate_random_q_element(height, n);
    let d2 = generate_random_q_element(height, n);
    let d3 = generate_random_q_element(height, n);
    let b = generate_random_q_element(height, n);
    let a = generate_random_q_element(height, n);
    let w_i = compute_w_i();
    SetupParams {
        constraints,
        height,
        n,
        q, q64,
        plans,
        g_matrix,
        d0, d1, d2, d3, b, a,
        w_i,
    }
}