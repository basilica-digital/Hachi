mod math;
mod mlp;
mod sumcheck;
mod utils;
mod prep;
mod hachi;

use tfhe_ntt::prime32::Plan;

use std::slice;
use std::fs::File;
use std::ptr::NonNull;
use std::time::Instant;
use std::arch::x86_64::*;

use rayon::prelude::*;

use crate::math::field_simd::*;
use crate::math::fields::*;
use crate::math::rns::*;
use crate::math::eq::*;
use crate::sumcheck::prove::sumcheck_prove;
use crate::utils::random::*;
use crate::utils::ds::*;
use crate::mlp::mle::{mle_w, mle_m,WBundle, MBundle};
use crate::hachi::commit::*;
use crate::prep::foldwitness::*;
use crate::hachi::setup::*;
use crate::utils::fs::*;
use crate::hachi::verify::*;
use crate::hachi::prove::*;

const Q: u32 = 4294967197;


#[allow(non_snake_case)]
fn main(){
    let start = Instant::now();
    
    // Setup
    let params = setup();
    
    let duration = start.elapsed();
    println!("Setup: {:?}", duration);
    let start = Instant::now();

    // prepare polynomial
    let s  = generate_random_data_4bit_packed(params.height*params.n, params.n);
    
    let duration = start.elapsed();
    println!("Prepare polynomial f: {:?}", duration);
    let start = Instant::now();
    
    // commit
    let commitment = Commit(&params, &s);

    let duration = start.elapsed();
    println!("Commit: {:?}", duration);
    let start = Instant::now();

    // prove
    let (proof, target_sum, tau_0, alpha) = unsafe { prove(&params, &s, commitment) };
    
    let duration = start.elapsed();
    println!("Prove: {:?}", duration);
    let start_verify = std::time::Instant::now();
    
    // Verifier
    let is_valid = verify(&params, &proof, target_sum, &tau_0, &alpha);

    let verify_duration = start_verify.elapsed();
    assert!(is_valid, "Failed");
    println!("Verify: {:?}", verify_duration);
}