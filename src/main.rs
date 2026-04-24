mod math;
mod mlp;
mod sumcheck;
mod utils;
mod prep;
mod hachi;

use std::path::PathBuf;
use std::process;
use std::time::Instant;

use clap::{Parser, Subcommand};

use crate::hachi::commit::Commit;
use crate::hachi::prove::prove;
use crate::hachi::setup::setup_with_seed;
use crate::hachi::verify::{sample_challenge, verify};
use crate::utils::ds::*;
use crate::utils::random::generate_random_data_4bit_packed;
use crate::utils::serialize::{
    read_commitment, read_proof, read_witness, write_commitment, write_proof, write_witness,
};

// Re-exports used by submodules via `crate::`.
pub use crate::hachi::setup::SetupParams;
pub use crate::math::field_simd::u642u32;
pub use crate::math::fields::{ext_sub, extmul};

const Q: u32 = 4294967197;

/// Hachi — lattice-based polynomial commitment scheme.
#[derive(Parser)]
#[command(name = "hachi")]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Generate a random 4-bit packed witness matching the setup dimensions.
    GenWitness {
        /// Output witness file path.
        #[arg(short, long, default_value = "witness.bin")]
        output: PathBuf,

        /// Deterministic seed for public parameter generation (controls
        /// matrix dimensions only — witness bytes are drawn from the OS RNG).
        #[arg(short, long, default_value_t = 0)]
        seed: u64,
    },

    /// Commit to a witness polynomial.
    Commit {
        /// Path to the raw 4-bit packed witness file.
        input: PathBuf,

        /// Output commitment file path.
        #[arg(short, long, default_value = "commitment.bin")]
        output: PathBuf,

        /// Deterministic seed for public parameter generation.
        #[arg(short, long, default_value_t = 0)]
        seed: u64,
    },

    /// Generate a sumcheck proof for a commitment.
    Prove {
        /// Path to the raw 4-bit packed witness file.
        #[arg(short, long)]
        witness: PathBuf,

        /// Path to the commitment file produced by `commit`.
        #[arg(short, long)]
        commitment: PathBuf,

        /// Output proof file path.
        #[arg(short, long, default_value = "proof.bin")]
        output: PathBuf,

        /// Deterministic seed for public parameter generation.
        #[arg(short, long, default_value_t = 0)]
        seed: u64,
    },

    /// Verify a sumcheck proof against a commitment.
    Verify {
        /// Path to the commitment file.
        #[arg(short, long)]
        commitment: PathBuf,

        /// Path to the proof file produced by `prove`.
        #[arg(short, long)]
        proof: PathBuf,

        /// Deterministic seed for public parameter generation.
        #[arg(short, long, default_value_t = 0)]
        seed: u64,
    },
}

#[allow(non_snake_case)]
fn main() {
    let cli = Cli::parse();

    match cli.command {
        Commands::GenWitness { output, seed } => {
            let start = Instant::now();
            let params = setup_with_seed(seed);
            println!("Setup: {:?}", start.elapsed());

            let start = Instant::now();
            // Witness shape matches what `prep::commit::commit` expects: a
            // `height_4 × n` grid of ring elements, each of `n` 4-bit
            // coefficients (i.e. `height_4 * n * n` nibbles total).
            let witness =
                generate_random_data_4bit_packed(params.height_4 * params.n, params.n);
            println!("Generate witness: {:?}", start.elapsed());

            write_witness(&output, &witness).unwrap_or_else(|e| {
                eprintln!("error: failed to write witness: {e}");
                process::exit(1);
            });
            println!("Witness written to {}", output.display());
        }

        Commands::Commit {
            input,
            output,
            seed,
        } => {
            let start = Instant::now();
            let params = setup_with_seed(seed);
            println!("Setup: {:?}", start.elapsed());

            let start = Instant::now();
            let s = read_witness(&input).unwrap_or_else(|e| {
                eprintln!("error: failed to read witness file: {e}");
                process::exit(1);
            });
            println!("Read witness: {:?}", start.elapsed());

            let start = Instant::now();
            let commitment = Commit(&params, &s);
            println!("Commit: {:?}", start.elapsed());

            write_commitment(&output, &commitment).unwrap_or_else(|e| {
                eprintln!("error: failed to write commitment: {e}");
                process::exit(1);
            });
            println!("Commitment written to {}", output.display());
        }

        Commands::Prove {
            witness,
            commitment,
            output,
            seed,
        } => {
            let start = Instant::now();
            let params = setup_with_seed(seed);
            println!("Setup: {:?}", start.elapsed());

            let start = Instant::now();
            let s = read_witness(&witness).unwrap_or_else(|e| {
                eprintln!("error: failed to read witness file: {e}");
                process::exit(1);
            });
            println!("Read witness: {:?}", start.elapsed());

            let start = Instant::now();
            let comm = read_commitment(&commitment).unwrap_or_else(|e| {
                eprintln!("error: failed to read commitment: {e}");
                process::exit(1);
            });
            println!("Read commitment: {:?}", start.elapsed());

            let challenge = sample_challenge(&params, seed);

            let start = Instant::now();
            let (proof, target_sum) = unsafe { prove(&params, &s, comm, &challenge) };
            println!("Prove: {:?}", start.elapsed());

            write_proof(&output, &proof, target_sum).unwrap_or_else(|e| {
                eprintln!("error: failed to write proof: {e}");
                process::exit(1);
            });
            println!("Proof written to {}", output.display());
        }

        Commands::Verify {
            commitment: _,
            proof,
            seed,
        } => {
            let start = Instant::now();
            let params = setup_with_seed(seed);
            println!("Setup: {:?}", start.elapsed());

            let start = Instant::now();
            let bundle = read_proof(&proof).unwrap_or_else(|e| {
                eprintln!("error: failed to read proof: {e}");
                process::exit(1);
            });
            println!("Read proof: {:?}", start.elapsed());

            let challenge = sample_challenge(&params, seed);

            let start = Instant::now();
            let is_valid = verify(
                &params,
                &bundle.proof,
                bundle.target_sum,
                &challenge,
            );
            println!("Verify: {:?}", start.elapsed());

            if is_valid {
                println!("VALID");
            } else {
                eprintln!("INVALID");
                process::exit(1);
            }
        }
    }
}