flags := "-C target-cpu=native -C codegen-units=1"
threads := "16"
seed := "0"
witness := "witness.bin"
commitment := "commitment.bin"
proof := "proof.bin"

# Shared cargo invocation: release build with native CPU tuning + Rayon threads.
cargo := 'RAYON_NUM_THREADS=' + threads + ' RUSTFLAGS="' + flags + '" cargo run --release --'

# Default recipe: list available commands.
default:
    @just --list

# Full pipeline: generate a witness, commit, prove, and verify.
all: gen-witness commit prove verify

# Generate a random witness sized to the setup parameters.
gen-witness:
    {{cargo}} gen-witness --seed {{seed}} --output {{witness}}

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
