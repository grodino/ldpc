use std::time::Instant;

use itertools::Itertools;
use rand::{seq::SliceRandom, thread_rng};
use rayon::prelude::*;
use serde::Serialize;

use crate::channel::{BinaryErasureChannel, BinaryMemorylessChannel, Channel};
use crate::coding::{Code, GallagerGenerator, Generate, Ldpc, MacKayNealGenerator};
use crate::decoding::{BECDecoder, BMCDecoder, Decoder};
use crate::utils::{binom, count_differences};

/// Result of a Montecarlo simulation of a code with varying snr (represented by
/// the values of `noise_parameter`)
#[derive(Debug, Serialize)]
pub struct MontecarloResult<'c, C: Code<N>, const N: usize> {
    name: String,
    noise_parameter: Vec<f64>,
    bit_error_rate: Vec<f64>,
    block_error_rate: Vec<f64>,
    code_rate: f64,
    decoding_time: f64,
    #[serde(borrow)]
    code: &'c C,
}

/// Simulate a Gallager code family
pub fn gallager_family_montecarlo<const N: usize, const M: usize>(
    snr_values: &Vec<f64>,
    n_realizations: usize,
    n_degree: usize,
) {
    // Dirty way to generate all possible degree distributions
    let degrees = (1..N)
        .cartesian_product(1..M)
        .filter(|&(bit_degree, parity_degree)| M * parity_degree == N * bit_degree)
        .take(n_degree);

    print!("[");
    for (i, (bit_degree, parity_degree)) in degrees.enumerate() {
        if i > 0 {
            print!(",");
        }
        eprintln!("EXP d_v = {bit_degree}, d_c = {parity_degree}");

        let code_generator = GallagerGenerator::<N, M>::new(bit_degree, parity_degree).unwrap();
        eprint!("\tCODE GENERATION ...");
        let start = Instant::now();

        let ldpc = code_generator.generate().unwrap();
        let duration = start.elapsed();
        eprintln!(" DONE in {:.6} s", duration.as_secs_f64());

        eprint!("\tMONTECARLO SIMULATION ...");
        let start = Instant::now();

        let result = montecarlo::<BinaryMemorylessChannel, Ldpc<N, M>, BMCDecoder<N, M>, N, M>(
            // let result = montecarlo::<BinaryErasureChannel, Ldpc<N, M>, BECDecoder<N, M>, N, M>(
            format!("d_v = {bit_degree}, d_c = {parity_degree}"),
            &ldpc,
            &snr_values[..],
            n_realizations,
        );
        let duration = start.elapsed().as_secs_f64();
        eprintln!(
            " DONE in {:.6} s ({:.6} s per codeword)",
            duration,
            duration / n_realizations as f64
        );
        print!("{}", serde_json::to_string(&result).unwrap());
    }
    print!("]");
}

/// Simulate a MacKay & Neal code family
pub fn mackay_family_montecarlo<'a, const N: usize, const M: usize>(
    snr_values: &Vec<f64>,
    n_realizations: usize,
    n_degree: usize,
    regular: bool,
    full_rank: bool,
    n_perturbations: usize,
) {
    let generators =
        (1..N) // Possible degree values for check nodes
            .cartesian_product(1..M) // Possible degree values for bit nodes
            .filter(|&(bit_degree, parity_degree)| M * parity_degree == N * bit_degree)
            .map(|(avg_bit_degree, avg_parity_degree)| {
                let n_retry = 2 * binom(N, avg_bit_degree);

                match regular {
                    true => MacKayNealGenerator::<N, M>::new_regular(
                        avg_bit_degree,
                        avg_parity_degree,
                        2 * binom(N, avg_bit_degree),
                        1 + N / 3,
                    )
                    .unwrap(),

                    false => {
                        let mut parity_degree = [avg_parity_degree; M];
                        let mut bit_degree = [avg_bit_degree; N];

                        // Add a total of n_perturbations^2 links
                        bit_degree[..n_perturbations].fill(avg_bit_degree + 3);
                        // Have to add them to the parity bits to preserve the total
                        // degree conservation equality
                        parity_degree[..n_perturbations].fill(avg_parity_degree + 3);

                        // shuffle the degree distributions
                        let mut rng = thread_rng();
                        bit_degree.shuffle(&mut rng);
                        parity_degree.shuffle(&mut rng);

                        MacKayNealGenerator::<N, M>::new(
                            bit_degree,
                            parity_degree,
                            full_rank,
                            n_retry,
                            1_000,
                            1 + N / 3,
                        )
                        .unwrap()
                    }
                }
            })
            .take(n_degree);

    print!("[");
    for (i, code_generator) in generators.enumerate() {
        if i > 0 {
            print!(",");
        }
        eprintln!("EXP {}", code_generator);

        eprint!("\tCODE GENERATION ...");
        let start = Instant::now();
        let ldpc = code_generator.generate().unwrap();
        let duration = start.elapsed();
        eprintln!(" DONE in {:.6} s", duration.as_secs_f64());

        eprint!("\tMONTECARLO SIMULATION ...");
        let result = montecarlo::<BinaryMemorylessChannel, Ldpc<N, M>, BMCDecoder<N, M>, N, M>(
            format!("{}", code_generator),
            &ldpc,
            &snr_values[..],
            n_realizations,
        );
        eprintln!(
            " DONE in {:.6} s ({:.6} s per codeword)",
            result.decoding_time * (n_realizations * snr_values.len()) as f64,
            result.decoding_time,
        );
        print!("{}", serde_json::to_string(&result).unwrap());
    }
    print!("]");
}

/// Montecarlo simulation of a given (channel, code, decoder) system.
pub fn montecarlo<'c, CH, CO, DE, const N: usize, const M: usize>(
    name: String,
    code: &'c CO,
    snr_values: &[f64],
    n_samples: usize,
) -> MontecarloResult<'c, CO, N>
where
    CH: Channel + std::marker::Sync,
    CO: 'c + Code<N>,
    DE: Decoder<'c, CH, CO, N> + std::marker::Sync,
{
    let message = CO::all_zero();
    let sent = CH::all_zero();
    let start = Instant::now();

    let (bit_error_rate, block_error_rate): (Vec<_>, Vec<_>) = snr_values
        .iter()
        .map(|&snr| {
            let channel: CH = CH::from_snr(snr);
            let decoder = DE::new(&channel, &code);

            // Compute bit error rate and block error rate
            let (bit_error, block_error) = (0..n_samples)
                .into_par_iter() // Run simulations in parallel
                .map(|_| {
                    let mut rng = thread_rng();

                    let received = channel.transmit(&sent, &mut rng);
                    let decoded = decoder.decode(&received);
                    let n_differences = count_differences(&message, &decoded);

                    (
                        n_differences as f64 / N as f64,
                        if n_differences > 0 { 1u64 } else { 0 },
                    )
                })
                .reduce(
                    || (0.0, 0u64),
                    |(acc_bit, acc_block), (bit_error, block_error)| {
                        (acc_bit + bit_error, acc_block + block_error)
                    },
                );

            (
                bit_error / n_samples as f64,
                block_error as f64 / n_samples as f64,
            )
        })
        .unzip();

    MontecarloResult {
        name: name,
        noise_parameter: snr_values.into(),
        bit_error_rate,
        block_error_rate,
        code_rate: (N - M) as f64 / N as f64,
        code,
        decoding_time: start.elapsed().as_secs_f64() / (n_samples * snr_values.len()) as f64,
    }
}
