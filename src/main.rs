//! Generate LDPC codes, encode code words and decode with message passing
use clap::{Parser, Subcommand};

mod channel;
mod coding;
mod decoding;
mod experiments;
mod utils;

use experiments::*;

/// Run experiments with LDPC codes
///
/// The output is a JSON object to ease the data processing. Each experiment is
/// identified by its name. Each contain a list code sizes. For each code size,
/// we generate code instances and simulate them `n_realizations` times. The
/// averaged results are then returned.
#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// (start, end, n) : the SNR ranges from `start` dB to `end` dB (inclusive)
    /// with n log spaced points. Example: --snr="-1,2,10"
    #[clap(long, parse(try_from_str = parse_tuple))]
    snr: Option<(f64, f64, usize)>,

    /// Number of channel realizations to average the values on
    #[clap(long, short = 'r', default_value_t = 1_000)]
    n_realizations: usize,

    /// Number of degree distributions to consider
    #[clap(long, short = 'd', default_value_t = 1)]
    n_degree: usize,

    /// The experiment to run
    #[clap(subcommand)]
    experiment: Experiment,
}

// #[derive(Debug, ArgEnum, Clone, PartialEq, Eq, PartialOrd, Ord)]
#[derive(Subcommand, Debug)]
enum Experiment {
    /// Generate a set of Gallager codes with same length and size but different
    /// degrees and return their averaged performance on BEC and BMC channel
    Gallager,

    /// Generate a set of regular MacKay & Neal codes with same length and size
    /// but different degrees and return their averaged performance on BEC and
    /// BMC channel
    MacKay {
        /// If flag set, generate only regular degree distributions
        #[clap(short, long, parse(from_flag))]
        regular: bool,

        /// If flag set, generate full rank parity check matrices (only
        /// applicable to non-regular degree distributions here)
        #[clap(short, long, parse(from_flag))]
        full_rank: bool,
    },
}

fn parse_tuple(
    s: &str,
) -> Result<(f64, f64, usize), Box<dyn std::error::Error + Send + Sync + 'static>> {
    let mut chunks = s.split(',');
    Ok((
        chunks
            .next()
            .ok_or_else(|| "error parsing tuple")?
            .parse()?,
        chunks
            .next()
            .ok_or_else(|| "error parsing tuple")?
            .parse()?,
        chunks
            .next()
            .ok_or_else(|| "error parsing tuple")?
            .parse()?,
    ))
}

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let args = Args::parse();

    let snr_values = match args.snr {
        Some((start, end, n)) => (0..n)
            .map(|i| 10f64.powf(start / 10.0 + i as f64 * (end - start) / (10.0 * n as f64)))
            .collect(),
        None => vec![0.1, 1.12883789, 2.97635144, 4.83293024, 7.8475997, 10.0],
    };

    println!("{{");

    match args.experiment {
        ////////////////////////////////////////////////////////////////////////////
        // EXPERIMENT 1: Gallager codes                                           //
        ////////////////////////////////////////////////////////////////////////////
        Experiment::Gallager => {
            println!(r#""name": "Gallager", "#);
            println!(r#""codes": ["#);

            // Rate 1/2 code
            gallager_family_montecarlo::<768, 384>(&snr_values, args.n_realizations, args.n_degree);
            println!(",");

            // Rate 2/3 code
            gallager_family_montecarlo::<768, 256>(&snr_values, args.n_realizations, args.n_degree);

            println!("\n]");
        }

        ////////////////////////////////////////////////////////////////////////////
        // EXPERIMENT 2: MacKay&Neal codes                                        //
        ////////////////////////////////////////////////////////////////////////////
        Experiment::MacKay { regular, full_rank } => {
            if regular && full_rank {
                unimplemented!("Cannot generate full rank, regular LDPC parity matrices for now");
            }
            println!(r#""name": "MacKay", "#);
            println!(r#""codes": ["#);

            // Rate 1/2 code
            mackay_family_montecarlo::<192, 96>(
                &snr_values,
                args.n_realizations,
                args.n_degree,
                regular,
                full_rank,
                4,
            );
            // mackay_family_montecarlo::<12, 6>(
            //     &snr_values,
            //     args.n_realizations,
            //     args.n_degree,
            //     regular,
            //     full_rank,
            //     2,
            // );

            println!(",");
            // mackay_family_montecarlo::<384, 192>(
            //     &snr_values,
            //     args.n_realizations,
            //     args.n_degree,
            //     regular,
            //     full_rank,
            // );
            // println!(",");
            // mackay_family_montecarlo::<768, 384>(
            //     &snr_values,
            //     args.n_realizations,
            //     args.n_degree,
            //     regular,
            //     full_rank,
            // );
            // println!(",");

            // Rate 2/3 codes
            mackay_family_montecarlo::<192, 64>(
                &snr_values,
                args.n_realizations,
                args.n_degree,
                regular,
                full_rank,
                192 / 3,
            );
            // println!(",");
            // mackay_family_montecarlo::<384, 256>(
            //     &snr_values,
            //     args.n_realizations,
            //     args.n_degree,
            //     regular,
            //     full_rank,
            // );
            // println!(",");
            // mackay_family_montecarlo::<768, 512>(
            //     &snr_values,
            //     args.n_realizations,
            //     args.n_degree,
            //     regular,
            //     full_rank,
            // );

            println!("\n]");
        }
    }

    println!("}}");

    Ok(())
}
