use rand_chacha::ChaCha12Rng;
use rand::SeedableRng;

pub struct FastTranscript {
    hasher: blake3::Hasher,
}

impl FastTranscript {
    pub fn new(label: &[u8]) -> Self {
        let mut hasher = blake3::Hasher::new();
        hasher.update(label);
        Self { hasher }
    }

    pub fn append_u32_slice(&mut self, label: &[u8], data: &[u32]) {
        self.hasher.update(label);
        let bytes: &[u8] = bytemuck::cast_slice(data);
        self.hasher.update(bytes);
    }

    pub fn get_challenge_rng(&self) -> ChaCha12Rng {
        let hash = self.hasher.finalize();
        let mut seed = [0u8; 32];
        seed.copy_from_slice(hash.as_bytes());
        ChaCha12Rng::from_seed(seed)
    }
}