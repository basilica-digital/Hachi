use std::path::PathBuf;
use std::process;
use std::time::Instant;

use clap::{Parser, Subcommand};

use hachi::hachi::commit::Commit;
use hachi::hachi::prove::prove;
use hachi::hachi::setup::setup_with_seed;
use hachi::hachi::verify::{sample_challenge, verify};
use hachi::utils::random::{generate_random_bytes_4bit_packed, stream_random_witness_to};
use hachi::utils::serialize::{
    read_commitment, read_proof, read_witness, write_commitment, write_proof, write_witness,
};
use hachi::utils::size::{format_size, parse_size};

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

        /// Target witness size, e.g. `1GB`, `10GiB`, `1TB`, `10KB`, or `1024`.
        /// Units use powers of 1024 and are case-insensitive. Must be a
        /// multiple of 16 bytes. When omitted, defaults to the size expected
        /// by the configured scheme parameters.
        ///
        /// NOTE: `commit`, `prove`, and `verify` still assume the
        /// scheme-default witness size; `--size` is useful for isolating
        /// `gen-witness` / witness-I/O benchmarks at other sizes.
        #[arg(short = 'S', long)]
        size: Option<String>,
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
        Commands::GenWitness {
            output,
            seed,
            size,
        } => {
            let start = Instant::now();
            let params = setup_with_seed(seed);
            println!("Setup: {:?}", start.elapsed());

            // Witness shape expected by `prep::commit::commit`: a `height_4 × n`
            // grid of ring elements, each with `n` 4-bit coefficients.
            let expected_bytes = params.height_4 * params.n * params.n / 2;

            let requested_bytes: usize = match size {
                Some(s) => parse_size(&s)
                    .and_then(|b| {
                        usize::try_from(b).map_err(|_| format!("size {b} exceeds usize"))
                    })
                    .unwrap_or_else(|e| {
                        eprintln!("error: invalid --size: {e}");
                        process::exit(1);
                    }),
                None => expected_bytes,
            };

            if requested_bytes % 16 != 0 {
                eprintln!(
                    "error: witness size must be a multiple of 16 bytes (got {requested_bytes})"
                );
                process::exit(1);
            }

            if requested_bytes != expected_bytes {
                eprintln!(
                    "warning: requested witness size {} differs from the scheme-default {}. \
                     `commit`, `prove`, and `verify` assume the scheme-default size and will \
                     fail or produce invalid results on a differently sized witness.",
                    format_size(requested_bytes as u64),
                    format_size(expected_bytes as u64),
                );
            }

            // Stream to disk when the witness would be large enough that
            // holding it in RAM is risky. The in-memory path keeps driving
            // the scheme-default (~4 GiB) case so tests/benchmarks are
            // unchanged.
            const STREAM_THRESHOLD: usize = 16 * (1 << 30); // 16 GiB
            const STREAM_CHUNK: usize = 256 * (1 << 20); // 256 MiB per chunk

            let start = Instant::now();
            if requested_bytes >= STREAM_THRESHOLD {
                stream_random_witness_to(&output, requested_bytes as u64, STREAM_CHUNK)
                    .unwrap_or_else(|e| {
                        eprintln!("error: failed to write witness: {e}");
                        process::exit(1);
                    });
                println!(
                    "Generate+write witness ({}, streamed): {:?}",
                    format_size(requested_bytes as u64),
                    start.elapsed()
                );
            } else {
                let witness = generate_random_bytes_4bit_packed(requested_bytes);
                println!(
                    "Generate witness ({}): {:?}",
                    format_size(requested_bytes as u64),
                    start.elapsed()
                );
                write_witness(&output, &witness).unwrap_or_else(|e| {
                    eprintln!("error: failed to write witness: {e}");
                    process::exit(1);
                });
            }
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