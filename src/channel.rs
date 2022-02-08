//! Physical channel description and simulation
use rand::Rng;
use rand_distr::{Bernoulli, Distribution, Normal};

/// A stationnary channel, determined only by its SNR.
pub trait Channel {
    /// The input type of the channel (usually binary)
    type Input: Copy + PartialOrd + Sync + Default;
    /// The output type of the channel (usally binary, binary + erasure or continuous)
    type Output: Copy + PartialOrd + Sync + Default;

    /// How the all zero codeword would be represented as the input of this channel
    fn all_zero<const N: usize>() -> [Self::Input; N];

    /// Transmit a sequence of N bits through the channel
    fn transmit<R: Rng, const N: usize>(
        &self,
        sequence: &[Self::Input; N],
        rng: &mut R,
    ) -> [Self::Output; N];

    /// Instanciate the Channel with the given SNR per bit assuming the signal
    /// has a unitary energy per bit.
    fn from_snr(snr: f64) -> Self;
}

pub struct BinaryErasureChannel {
    pub erasure_proba: f64,
}

impl Channel for BinaryErasureChannel {
    type Input = u8;
    type Output = u8;

    fn transmit<R: Rng, const N: usize>(
        &self,
        sequence: &[Self::Input; N],
        rng: &mut R,
    ) -> [Self::Output; N] {
        let bernoulli = Bernoulli::new(self.erasure_proba).unwrap();

        let mut received = [0u8; N];
        received.copy_from_slice(
            &bernoulli
                .sample_iter(rng)
                .take(N)
                .zip(sequence)
                .map(|(error, &symbol)| match error {
                    true => 2,
                    false => symbol,
                })
                .collect::<Vec<u8>>(),
        );

        received
    }

    fn from_snr(snr: f64) -> Self {
        Self {
            erasure_proba: (1.0 / snr).min(1.0),
        }
    }

    fn all_zero<const N: usize>() -> [Self::Input; N] {
        [0u8; N]
    }
}

/// A Gaussian binary memoryless channel with given noise variance.
pub struct BinaryMemorylessChannel {
    pub noise_variance: f64,
}

impl Channel for BinaryMemorylessChannel {
    type Input = f64;
    type Output = f64;

    fn transmit<R: Rng, const N: usize>(
        &self,
        sequence: &[Self::Input; N],
        rng: &mut R,
    ) -> [Self::Output; N] {
        let gaussian = Normal::new(0.0, self.noise_variance).unwrap();

        let mut received = [0f64; N];
        received.copy_from_slice(
            &gaussian
                .sample_iter(rng)
                .take(N)
                .zip(sequence)
                .map(|(noise, &symbol)| symbol + noise)
                .collect::<Vec<f64>>(),
        );

        received
    }

    /// Instanciate a `BinaryMemorylessChannel` with the given SNR per bit
    /// assuming the signal has a unitary energy per bit.
    fn from_snr(snr: f64) -> Self {
        Self {
            noise_variance: 1.0 / snr,
        }
    }

    fn all_zero<const N: usize>() -> [Self::Input; N] {
        [1f64; N]
    }
}
