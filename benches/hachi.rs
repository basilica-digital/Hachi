//! Criterion benchmarks for the Hachi commitment / prove / verify pipeline.
//!
//! Benchmarks are parameterised over witness size in bytes. Today the scheme
//! only supports the default 4 GiB witness — adding a new size to
//! `WITNESS_SIZES` will fail at setup until the inner SIMD loops are
//! generalised over `height_2` / `height_4`. The constant lives at the top
//! of the file so it is the single place to extend once larger sizes work.
//!
//! Run with:
//! ```sh
//! RUSTFLAGS="-C target-cpu=native" cargo bench --bench hachi
//! ```
//! Each fixture allocates ~4 GiB for the witness plus several GiB of
//! intermediate state (`prove` is the heaviest). Make sure the host has
//! enough RAM and disk-backed swap before invoking. The fixture is built
//! once per witness size and shared across `commit`, `prove`, and `verify`.

use std::collections::HashMap;
use std::sync::{Mutex, OnceLock};
use std::time::Duration;

use criterion::{
    criterion_group, criterion_main, BenchmarkId, Criterion, Throughput,
};

use hachi::hachi::commit::{Commit, Commitment};
use hachi::hachi::prove::{prove, SumcheckProof};
use hachi::hachi::setup::{setup_with_seed_and_size, SetupParams, DEFAULT_WITNESS_SIZE_BYTES};
use hachi::hachi::verify::{sample_challenge, verify, Challenge};
use hachi::utils::ds::AlignedU8Vec;
use hachi::utils::random::generate_random_bytes_4bit_packed;
use hachi::utils::size::format_size;

/// Witness sizes (in bytes) to benchmark. Currently only the scheme-default
/// 4 GiB is wired up; extend this list as the underlying setup learns to
/// scale `height_2` / `height_4` to other sizes.
const WITNESS_SIZES: &[u64] = &[DEFAULT_WITNESS_SIZE_BYTES];

const SEED: u64 = 0;

/// One-shot fixture: setup parameters + witness + Fiat-Shamir challenge.
/// Reused across every commit/prove/verify benchmark sample for a given
/// witness size — building it costs minutes and ~4 GiB so we want one copy.
struct Fixture {
    params: SetupParams,
    witness: AlignedU8Vec,
    challenge: Challenge,
    /// Pre-computed proof material, populated lazily on first verify call so
    /// the commit/prove benches don't pay for it.
    verify_artifacts: OnceLock<(SumcheckProof, [u32; 4])>,
}

impl Fixture {
    fn new(size_bytes: u64) -> Self {
        let params = setup_with_seed_and_size(SEED, size_bytes)
            .expect("witness size unsupported by current setup");
        let expected = params.witness_size_bytes();
        assert_eq!(
            expected, size_bytes,
            "setup produced a witness size that doesn't match the requested one",
        );

        let witness_bytes: usize = expected
            .try_into()
            .expect("witness size exceeds usize on this target");
        let witness = generate_random_bytes_4bit_packed(witness_bytes);
        let challenge = sample_challenge(&params, SEED);

        Fixture {
            params,
            witness,
            challenge,
            verify_artifacts: OnceLock::new(),
        }
    }

    fn verify_artifacts(&self) -> &(SumcheckProof, [u32; 4]) {
        self.verify_artifacts.get_or_init(|| {
            let commitment = Commit(&self.params, &self.witness);
            unsafe { prove(&self.params, &self.witness, commitment, &self.challenge) }
        })
    }
}

/// Shared fixture cache keyed by witness size. Holding `&'static Fixture`
/// lets the bench closures borrow it for the lifetime of the program
/// without juggling lifetimes through criterion's API. Leaking is fine
/// here — the process exits after the bench run.
fn fixture(size_bytes: u64) -> &'static Fixture {
    static CACHE: OnceLock<Mutex<HashMap<u64, &'static Fixture>>> = OnceLock::new();
    let cache = CACHE.get_or_init(|| Mutex::new(HashMap::new()));
    let mut guard = cache.lock().expect("fixture cache poisoned");
    if let Some(&f) = guard.get(&size_bytes) {
        return f;
    }
    let leaked: &'static Fixture = Box::leak(Box::new(Fixture::new(size_bytes)));
    guard.insert(size_bytes, leaked);
    leaked
}

fn bench_commit(c: &mut Criterion) {
    let mut group = c.benchmark_group("commit");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(60));

    for &size in WITNESS_SIZES {
        let fix = fixture(size);
        group.throughput(Throughput::Bytes(size));
        group.bench_with_input(
            BenchmarkId::from_parameter(format_size(size)),
            &size,
            |b, _| {
                b.iter(|| Commit(&fix.params, &fix.witness));
            },
        );
    }

    group.finish();
}

/// `prove` consumes its `Commitment`, so we rebuild one per iteration via
/// `iter_batched` — the rebuild cost is excluded from the timed region.
fn bench_prove(c: &mut Criterion) {
    let mut group = c.benchmark_group("prove");
    group.sample_size(10);
    group.measurement_time(Duration::from_secs(120));

    for &size in WITNESS_SIZES {
        let fix = fixture(size);
        group.throughput(Throughput::Bytes(size));
        group.bench_with_input(
            BenchmarkId::from_parameter(format_size(size)),
            &size,
            |b, _| {
                b.iter_batched(
                    || -> Commitment { Commit(&fix.params, &fix.witness) },
                    |commitment| unsafe {
                        prove(&fix.params, &fix.witness, commitment, &fix.challenge)
                    },
                    criterion::BatchSize::PerIteration,
                );
            },
        );
    }

    group.finish();
}

/// Verify is cheap relative to commit/prove, so we use a larger sample size.
/// One valid proof is built lazily via `Fixture::verify_artifacts` and
/// re-verified for each sample.
fn bench_verify(c: &mut Criterion) {
    let mut group = c.benchmark_group("verify");
    group.sample_size(50);
    group.measurement_time(Duration::from_secs(30));

    for &size in WITNESS_SIZES {
        let fix = fixture(size);
        let (proof, target_sum) = fix.verify_artifacts();
        group.throughput(Throughput::Bytes(size));
        group.bench_with_input(
            BenchmarkId::from_parameter(format_size(size)),
            &size,
            |b, _| {
                b.iter(|| verify(&fix.params, proof, *target_sum, &fix.challenge));
            },
        );
    }

    group.finish();
}

criterion_group!(benches, bench_commit, bench_prove, bench_verify);
criterion_main!(benches);
