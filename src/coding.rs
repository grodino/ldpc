//! Error correcting code description and generation

use rand::{seq::index, thread_rng};
use serde::Serialize;
use std::fmt;

use crate::utils::{print_dense, rank};

/// The description of a (linear) code (length, parity check matrix ...)
///
/// N is the code length (the length of the codewords)
pub trait Code<const N: usize> {
    /// The output type of the codewords
    type Output: Copy + PartialOrd + Sync;

    /// The all zero codeword
    fn all_zero() -> [Self::Output; N];
}

/// A code generator allows to generate the given code with a given structure
pub trait Generate<C: Code<N>, const N: usize> {
    fn generate(&self) -> Result<C, GeneratorError>;
}

/// Binary LDPC code
#[derive(Debug, Serialize)]
pub struct Ldpc<const N: usize, const M: usize> {
    /// Length of the code (number of symbols in a codeword)
    // length: N,
    /// Size of the codebook (number of different codewords)
    // size: usize,
    /// Dimension of the code (rank of the encoder)
    pub dimension: usize,
    /// LDPC parity check matrix $H$ stored in Binary Sparse Row format
    /// ```parity_row[bit_node][j] = check_node``` TODO: remove pub
    pub parity_row: Vec<Vec<usize>>,
    /// LDPC parity check matrix $H$ stored in Binary Sparse Column format TODO: remove pub
    pub parity_col: Vec<Vec<usize>>,
    /// The total degree of the Tanner graph (ie the number of ones in the
    /// parity matrix H). It is equal to the number of edge out of all the bit
    /// nodes (which is equal to the number of edges out of all the check nodes)
    pub total_degree: u64,
    #[serde(skip)]
    pub parity_check: [[u8; N]; M],
}

impl<const N: usize, const M: usize> Ldpc<N, M> {
    /// Instanciate an LDPC code from its dense parity check matrix
    pub fn from_dense_parity_check(parity_check: [[u8; N]; M]) -> Self {
        // Store the matrix as a list of rows represented by non-zero
        // coefficients indices
        let parity_row = parity_check
            .iter()
            .map(|line| {
                // Keep the indexes of the non zero elements
                line.iter()
                    .enumerate()
                    .filter_map(|(i, &coeff)| {
                        if coeff > 0 {
                            return Some(i);
                        }

                        None
                    })
                    .collect()
            })
            .collect();

        // Store the matrix as a list of columns represented by their non-zero
        // elements
        let mut parity_col = Vec::with_capacity(M);

        for j in 0..N {
            let mut non_zero = Vec::with_capacity(N);

            for i in 0..M {
                if parity_check[i][j] != 0 {
                    // The check to bit are intially all erasures
                    non_zero.push(i);
                }
            }

            non_zero.shrink_to_fit();
            parity_col.push(non_zero);
        }

        Self {
            parity_check,
            dimension: rank(&parity_check),
            parity_row,
            parity_col,
            total_degree: parity_check
                .into_iter()
                .flatten()
                .fold(0u64, |acc, item| acc + item as u64) as u64,
        }
    }

    /// Returns the length N of the code
    pub const fn length(&self) -> usize {
        N
    }
}

impl<const N: usize, const M: usize> Code<N> for Ldpc<N, M> {
    type Output = u8;

    fn all_zero() -> [Self::Output; N] {
        [0u8; N]
    }
}

/// Generate a (d_v, d_c)-regular LDPC parity check matrix with gallager
/// generation method.
///
/// d_v is the bit degree and d_c is the parity check equation degree.
pub struct GallagerGenerator<const N: usize, const M: usize> {
    bit_degree: usize,
    parity_degree: usize,
}

impl<const N: usize, const M: usize> fmt::Display for GallagerGenerator<N, M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "Gallager, d_v = {}, d_c = {}",
            self.bit_degree, self.parity_degree
        )
    }
}

impl<const N: usize, const M: usize> GallagerGenerator<N, M> {
    pub fn new(bit_degree: usize, parity_degree: usize) -> Result<Self, GeneratorError> {
        if M * parity_degree != N * bit_degree || bit_degree == 0 {
            return Err(GeneratorError::DegreeDistributionError);
        }

        Ok(GallagerGenerator {
            bit_degree,
            parity_degree,
        })
    }
}

impl<const N: usize, const M: usize> Generate<Ldpc<N, M>, N> for GallagerGenerator<N, M> {
    fn generate(&self) -> Result<Ldpc<N, M>, GeneratorError> {
        let mut rng = thread_rng();
        let mut parity_check = [[0u8; N]; M];

        // Generate the first submatrix H1 of dimensions m/d_v × n such that the
        // i-th row, i ≤ m/d_v , has non-zero entries in the ((i − 1)d_c + 1)-th
        // to id_c-th columns
        for i in 0..(M / self.bit_degree) {
            for j in i * self.parity_degree..(i + 1) * self.parity_degree {
                parity_check[i][j] = 1;
            }
        }

        // Construct the rest of the lines by column permutation
        for i_submatrix in 1..self.bit_degree {
            let permutation = index::sample(&mut rng, N, N).into_vec();

            for i in 0..(M / self.bit_degree) {
                for j in 0..N {
                    parity_check[i_submatrix * (M / self.bit_degree) + i][j] =
                        parity_check[i][permutation[j]];
                }
            }
        }

        Ok(Ldpc::from_dense_parity_check(parity_check))
    }
}

/// Represents a node degree distribution, either regular or irregular
#[derive(Debug)]
enum DegreeDistribution<const K: usize> {
    Regular(usize),
    Distribution([usize; K]),
}

/// Generate a (d_v, d_c)-regular LDPC parity check matrix with MacKay & Neal
/// generation method.
///
/// d_v is the bit degree and d_c is the parity check equation degree. The
/// output code can be regular, irregular and forced to have full-rank parity
/// check
pub struct MacKayNealGenerator<const N: usize, const M: usize> {
    bit_degree: DegreeDistribution<N>,
    parity_degree: DegreeDistribution<M>,
    full_rank: Option<bool>,
    /// Number of retries for the generation of one column
    n_retries: usize,
    backtracking_depth: usize,
    n_rank_retries: usize,
}

impl<const N: usize, const M: usize> fmt::Display for MacKayNealGenerator<N, M> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let full_rank = self.full_rank.unwrap_or(false);

        match (&self.bit_degree, &self.parity_degree) {
            (
                DegreeDistribution::Regular(bit_degree),
                DegreeDistribution::Regular(parity_degree),
            ) => write!(f, "d_v = {}, d_c = {}", bit_degree, parity_degree),
            (
                DegreeDistribution::Distribution(bit_degree),
                DegreeDistribution::Distribution(parity_degree),
            ) => write!(
                f,
                "d_v = {:0.2}, d_c = {:0.2}",
                bit_degree.iter().sum::<usize>() as f64 / N as f64,
                parity_degree.iter().sum::<usize>() as f64 / M as f64,
            ),
            _ => unreachable!(),
        }
    }
}

impl<const N: usize, const M: usize> MacKayNealGenerator<N, M> {
    /// Instanciate a MacKay & Neal irregular degree LDPC code generator
    pub fn new(
        bit_degree: [usize; N],
        parity_degree: [usize; M],
        full_rank: bool,
        n_retries: usize,
        n_rank_retries: usize,
        backtracking_depth: usize,
    ) -> Result<Self, GeneratorError> {
        if parity_degree.iter().sum::<usize>() != bit_degree.iter().sum::<usize>() {
            return Err(GeneratorError::DegreeDistributionError);
        }

        Ok(Self {
            bit_degree: DegreeDistribution::Distribution(bit_degree),
            parity_degree: DegreeDistribution::Distribution(parity_degree),
            full_rank: Some(full_rank),
            n_retries,
            n_rank_retries,
            backtracking_depth,
        })
    }

    /// Instanciate a regular MacKay & Neal LDPC code generator
    pub fn new_regular(
        bit_degree: usize,
        parity_degree: usize,
        n_retries: usize,
        backtracking_depth: usize,
    ) -> Result<Self, GeneratorError> {
        if M * parity_degree != N * bit_degree {
            return Err(GeneratorError::DegreeDistributionError);
        }

        Ok(Self {
            bit_degree: DegreeDistribution::Regular(bit_degree),
            parity_degree: DegreeDistribution::Regular(parity_degree),
            full_rank: None,
            n_retries,
            n_rank_retries: 0,
            backtracking_depth,
        })
    }

    /// Generate a (d_v, d_c)-regular LDPC parity check matrix with MacKay and
    /// Neal's generation method. d_v is the bit degree and d_c is the parity
    /// check equation degree.
    fn generate_regular(
        &self,
        bit_degree: usize,
        parity_degree: usize,
    ) -> Result<Ldpc<N, M>, GeneratorError> {
        let mut rng = thread_rng();

        // The backtracking depth
        let t = self.backtracking_depth;
        // The max number of retries for the generation of one column
        let mut parity_check = [[0u8; N]; M];

        // Backtrack i.e erase the last t columns
        let backtrack = |parity_check: &mut [[u8; N]; M], current_col: usize| -> usize {
            let s = if current_col >= t { current_col - t } else { 0 };

            // erase the last t + 1 columns
            for l in s..current_col {
                for k in 0..M {
                    parity_check[k][l] = 0;
                }
            }

            s
        };

        // Check if the row degree profile is met until the given column
        let check_row_degree = |parity_check: &[[u8; N]; M], until: usize| -> bool {
            parity_check
                .iter()
                .map(|row| {
                    row.into_iter()
                        .take(until + 1)
                        .fold(0usize, |acc, &x| acc + x as usize)
                })
                .all(|degree| degree <= parity_degree)
        };

        // The index of the column being generated
        let mut j = 0;

        while j < N {
            let mut successful_generation = false;

            for _ in 0..self.n_retries {
                // Delete the column j
                for i in 0..M {
                    parity_check[i][j] = 0;
                }

                // Randomly generate column j with weight d_v
                let indices = index::sample(&mut rng, M, bit_degree).into_vec();
                for i in indices {
                    parity_check[i][j] = 1;
                }

                // Check the degree constraint
                if check_row_degree(&parity_check, j) {
                    successful_generation = true;
                    break;
                }
            }

            j = match successful_generation {
                true => j + 1,
                false => backtrack(&mut parity_check, j),
            }
        }

        Ok(Ldpc::from_dense_parity_check(parity_check))
    }

    /// Generate a (d_v, d_c)-irregular LDPC parity check matrix with MacKay and
    /// Neal's generation method. d_v is the bit degree distribution and d_c is
    /// the parity check equation degree distribution.
    fn generate_irregular(
        &self,
        bit_degree: &[usize; N],
        parity_degree: &[usize; M],
    ) -> Result<Ldpc<N, M>, GeneratorError> {
        let mut rng = thread_rng();

        // The backtracking depth
        let t = self.backtracking_depth;
        // The max number of retries for the generation of one column
        let mut parity_check = [[0u8; N]; M];

        // Backtrack i.e erase the last t columns
        let backtrack = |parity_check: &mut [[u8; N]; M], current_col: usize| -> usize {
            let s = if current_col >= t { current_col - t } else { 0 };

            // erase the last t + 1 columns
            for l in s..current_col {
                for k in 0..M {
                    parity_check[k][l] = 0;
                }
            }

            s
        };

        // Check if the row degree profile is met until the given column
        let check_row_degree = |parity_check: &[[u8; N]; M], until: usize| -> bool {
            parity_check
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    (
                        i,
                        row.into_iter()
                            .take(until + 1)
                            .fold(0usize, |acc, &x| acc + x as usize),
                    )
                })
                .all(|(i, degree)| degree <= parity_degree[i])
        };

        // The index of the column being generated
        let mut j = 0;

        while j < N {
            let mut successful_generation = false;

            for _ in 0..self.n_retries {
                // Delete the column j
                for i in 0..M {
                    parity_check[i][j] = 0;
                }

                // Randomly generate column j with weight d_v
                let indices = index::sample(&mut rng, M, bit_degree[j]).into_vec();
                for i in indices {
                    parity_check[i][j] = 1;
                }

                // Check the degree constraint
                if check_row_degree(&parity_check, j) {
                    successful_generation = true;
                    break;
                }
            }

            j = match successful_generation {
                true => j + 1,
                false => backtrack(&mut parity_check, j),
            }
        }

        // eprintln!();
        // print_dense(&parity_check);

        Ok(Ldpc::from_dense_parity_check(parity_check))
    }

    /// Generate a (d_v, d_c)-irregular LDPC parity check matrix with MacKay and
    /// Neal's generation method. d_v is the bit degree distribution and d_c is
    /// the parity check equation degree distribution.
    fn generate_irregular_full_rank(
        &self,
        bit_degree: &[usize; N],
        parity_degree: &[usize; M],
    ) -> Result<Ldpc<N, M>, GeneratorError> {
        let mut rng = thread_rng();

        // The backtracking depth
        let t = self.backtracking_depth;

        // Backtrack i.e erase the last t columns
        let backtrack = |parity_check: &mut [[u8; N]; M], current_col: usize| -> usize {
            let s = if current_col >= t { current_col - t } else { 0 };

            // erase the last t + 1 columns
            for l in s..current_col {
                for k in 0..M {
                    parity_check[k][l] = 0;
                }
            }

            s
        };

        // Check if the row degree profile is met until the given column
        let check_row_degree = |parity_check: &[[u8; N]; M], until: usize| -> bool {
            parity_check
                .iter()
                .enumerate()
                .map(|(i, row)| {
                    (
                        i,
                        row.into_iter()
                            .take(until + 1)
                            .fold(0usize, |acc, &x| acc + x as usize),
                    )
                })
                .all(|(i, degree)| degree <= parity_degree[i])
        };

        for _ in 0..self.n_rank_retries {
            let mut parity_check = [[0u8; N]; M];
            // The index of the column being generated
            let mut j = 0;

            while j < N {
                let mut successful_generation = false;

                for _ in 0..self.n_retries {
                    // Delete the column j
                    for i in 0..M {
                        parity_check[i][j] = 0;
                    }

                    // Randomly generate column j with weight d_v
                    let indices = index::sample(&mut rng, M, bit_degree[j]).into_vec();
                    for i in indices {
                        parity_check[i][j] = 1;
                    }

                    // Check the degree constraint
                    if check_row_degree(&parity_check, j) {
                        successful_generation = true;
                        break;
                    }
                }

                // If the degree constraint could not be met, we need to backtrack
                j = match successful_generation {
                    true => j + 1,
                    false => backtrack(&mut parity_check, j),
                }
            }

            // Check the rank constraint
            let r = rank(&parity_check);
            // println!("rank {r}");

            if r == M {
                return Ok(Ldpc::from_dense_parity_check(parity_check));
            }
        }

        Err(GeneratorError::NotFullRankError)
    }
}

impl<const N: usize, const M: usize> Generate<Ldpc<N, M>, N> for MacKayNealGenerator<N, M> {
    fn generate(&self) -> Result<Ldpc<N, M>, GeneratorError> {
        match (&self.bit_degree, &self.parity_degree) {
            (
                &DegreeDistribution::Regular(bit_degree),
                &DegreeDistribution::Regular(parity_degree),
            ) => self.generate_regular(bit_degree, parity_degree),
            (
                DegreeDistribution::Distribution(bit_degree),
                DegreeDistribution::Distribution(parity_degree),
            ) => match self.full_rank {
                Some(false) => self.generate_irregular(bit_degree, parity_degree),
                Some(true) => self.generate_irregular_full_rank(bit_degree, parity_degree),
                _ => unreachable!(),
            },
            _ => unreachable!(),
        }
    }
}

#[derive(Debug)]
/// A code generation error
pub enum GeneratorError {
    /// The degree distribution of the bit nodes and parity check nodes are not
    /// compatible
    DegreeDistributionError,
    /// When trying to generate the code, could not generate a full rank matrix
    NotFullRankError,
}

impl std::error::Error for GeneratorError {}

impl fmt::Display for GeneratorError {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        match self {
            GeneratorError::DegreeDistributionError => write!(f, "The degree distribution of the bit nodes and parity check nodes are not compatible"),
            GeneratorError::NotFullRankError => write!(f, "When trying to generate the code, could not generate a full rank matrix"),
        }
    }
}
