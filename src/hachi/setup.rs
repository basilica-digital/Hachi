use tfhe_ntt::prime32::Plan;
use crate::utils::random::{generate_random_q_element, generate_random_q_element_seeded};
use crate::math::fields::mod_inverse;
use crate::Q;
use crate::utils::ds::AlignedU32Vec;
use rand::SeedableRng;
use rand_chacha::ChaCha12Rng;

/// Witness size produced by the current scheme parameters: `height_4 * n * n / 2`
/// bytes, with `height_4 = 2^13` and `n = 2^10`.
pub const DEFAULT_WITNESS_SIZE_BYTES: u64 = (1u64 << 13) * (1u64 << 10) * (1u64 << 10) / 2;

pub struct SetupParams {
    pub constraints: usize,
    pub height_2: usize,
    pub height_4: usize,
    pub n: usize,
    pub q: u32,
    pub q64: u64,
    pub plans: Vec<Plan>,
    pub g_matrix_2: Vec<u32>,
    pub g_matrix_4: Vec<u32>,
    pub d: AlignedU32Vec,
    pub e: AlignedU32Vec,
    pub w_i: [u32; 7],
}

impl SetupParams {
    /// Witness size in bytes implied by the current matrix dimensions.
    ///
    /// 4-bit-packed: `height_4 * n * n / 2` bytes. Today this is 4 GiB; the
    /// constant is exposed so benchmarks and CLI callers can compute the
    /// supported size without duplicating the formula.
    pub fn witness_size_bytes(&self) -> u64 {
        (self.height_4 as u64) * (self.n as u64) * (self.n as u64) / 2
    }
}

fn compute_w_i() -> [u32; 7] {
    let mut w = [0u32; 7];
    for i in 0..7usize {
        let mut denom = 1u64;
        for j in 0..7usize {
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

/// Deterministic setup that targets a specific witness byte count.
///
/// Currently only [`DEFAULT_WITNESS_SIZE_BYTES`] (4 GiB) is supported because
/// many of the inner loops use literal `1<<13` / `1<<14` constants. The
/// signature is in place so callers (benchmarks, CLI) can sweep over
/// `witness_size_bytes` once the SIMD code below is generalised; for now
/// passing any other value returns an error.
pub fn setup_with_seed_and_size(
    seed: u64,
    witness_size_bytes: u64,
) -> Result<SetupParams, String> {
    if witness_size_bytes != DEFAULT_WITNESS_SIZE_BYTES {
        return Err(format!(
            "witness_size_bytes={witness_size_bytes} unsupported; only {} is wired up today",
            DEFAULT_WITNESS_SIZE_BYTES
        ));
    }
    Ok(setup_with_seed(seed))
}

/// Deterministic setup: same `seed` always produces identical parameters.
pub fn setup_with_seed(seed: u64) -> SetupParams {
    let constraints = 15;
    let height_2 = 1 << 14;
    let height_4 = 1 << 13;
    let n = 1 << 10;
    let q = 4294967197u32;
    let q64 = 4294967197u64;
    let moduli = [2079301633, 2079305729];
    let plans: Vec<Plan> = moduli.iter().map(|&m| Plan::try_new(1 << 11, m).unwrap()).collect();
    let mut g_matrix_2 = vec![0u32; 16];
    for i in 0..16 {
        g_matrix_2[i] = 1 << (i * 2);
    }
    let mut g_matrix_4 = vec![0u32; 8];
    for i in 0..8 {
        g_matrix_4[i] = 1 << (i * 4);
    }

    let mut master_rng = ChaCha12Rng::seed_from_u64(seed);
    let d = generate_random_q_element_seeded(height_2, n, &mut master_rng);
    let e = generate_random_q_element_seeded(64, n, &mut master_rng);

    let w_i = compute_w_i();
    SetupParams {
        constraints,
        height_2,
        height_4,
        n,
        q, q64,
        plans,
        g_matrix_2,
        g_matrix_4,
        d, e,
        w_i,
    }
}

pub fn setup() -> SetupParams{
    let constraints = 15;
    let height_2 = 1<<14;
    let height_4 = 1<<13;
    let n = 1<<10;
    let q = 4294967197u32;
    let q64 = 4294967197u64;
    let moduli = [2079301633, 2079305729];
    let plans: Vec<Plan> = moduli.iter().map(|&m| Plan::try_new(1<<11, m).unwrap()).collect();
    let mut g_matrix_2 = vec![0u32; 16];
    for i in 0..16{
        g_matrix_2[i] = 1<<(i*2);
    }
    let mut g_matrix_4 = vec![0u32; 8];
    for i in 0..8{
        g_matrix_4[i] = 1<<(i*4);
    }
    let d = generate_random_q_element(height_2, n);
    let e = generate_random_q_element(64, n);
    let w_i = compute_w_i();
    SetupParams{
        constraints,
        height_2,
        height_4,
        n,
        q, q64,
        plans,
        g_matrix_2,
        g_matrix_4,
        d, e,
        w_i,
    }
}