use crate::utils::ds::*;
use crate::prep::commit::*;
use crate::hachi::setup::SetupParams;

pub struct Commitment {
    pub u: AlignedU32Vec,
    pub r: AlignedU32Vec,
    pub t: AlignedU32Vec,
}

pub fn Commit(params: &SetupParams, s: &AlignedU8Vec) -> Commitment {
    let total_len = params.height_2 * params.n;
    let mut r = AlignedU32Vec {
        inner: vec![Align64([0u32; 16]); params.n / 16],
        len: params.n,
    };
    let mut t = AlignedU32Vec {
        inner: vec![Align64([0u32; 16]); total_len / 16],
        len: total_len,
    };
    let mut u = AlignedU32Vec {
        inner: vec![Align64([0u32; 16]); params.n / 16],
        len: params.n,
    };
    unsafe {
        commit(&mut u, &mut r[0..params.n], &mut t, s, &params.d);
    }
    Commitment{u, r, t}
}