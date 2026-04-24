flags := "-C target-cpu=native -C codegen-units=1"
threads := "16"
seed := "0"
# Override the witness byte size: e.g. `just size=1GB gen-witness`. Empty
# string means "use the size expected by the configured scheme parameters".
size := ""
witness := "witness.bin"
commitment := "commitment.bin"
proof := "proof.bin"

# Shared cargo invocation: release build with native CPU tuning + Rayon threads.
cargo := 'RAYON_NUM_THREADS=' + threads + ' RUSTFLAGS="' + flags + '" cargo run --release --'
size_arg := if size == "" { "" } else { "--size " + size }

# Default recipe: list available commands.
default:
    @just --list

# Full pipeline: generate a witness, commit, prove, and verify.
all: gen-witness commit prove verify

# Generate a random witness sized to the setup parameters (or `size` if set).
gen-witness:
    {{cargo}} gen-witness --seed {{seed}} --output {{witness}} {{size_arg}}

# Commit to the witness, producing a commitment file.
commit:
    {{cargo}} commit {{witness}} --seed {{seed}} --output {{commitment}}

# Generate a sumcheck proof for the commitment.
prove:
    {{cargo}} prove --witness {{witness}} --commitment {{commitment}} --seed {{seed}} --output {{proof}}

# Verify a sumcheck proof against the commitment.
verify:
    {{cargo}} verify --commitment {{commitment}} --proof {{proof}} --seed {{seed}}

# Remove generated witness, commitment, and proof artefacts.
clean:
    rm -f {{witness}} {{commitment}} {{proof}}

# Run criterion benchmarks (commit / prove / verify, parameterised by
# witness size). Each iteration allocates ~4 GiB plus the prover's
# intermediate state — make sure the host has the memory headroom.
bench:
    RAYON_NUM_THREADS={{threads}} RUSTFLAGS="{{flags}}" cargo bench --bench hachi
