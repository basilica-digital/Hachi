//! Hachi — lattice-based polynomial commitment scheme.
//!
//! Library entry point. Exposes the `commit`, `prove`, and `verify` pipeline
//! for downstream consumers (CLI binary, criterion benches, integration
//! tests). Internal modules keep their `crate::Q` / `crate::SetupParams` /
//! `crate::Align64*` references resolving here.

pub mod hachi;
pub mod math;
pub mod mlp;
pub mod prep;
pub mod sumcheck;
pub mod utils;

pub use crate::hachi::setup::SetupParams;
pub use crate::math::field_simd::u642u32;
pub use crate::math::fields::{ext_sub, extmul};
pub use crate::utils::ds::*;

pub const Q: u32 = 4294967197;
