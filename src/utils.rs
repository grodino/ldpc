//! Utilitary functions

/// Perform in place Gaussian elimination on a NxM matrix. It is assumed that
/// this matrix is full rank and that it its coefficients take binary values.
pub fn binary_elimination<const N: usize, const M: usize>(matrix: &mut [[u8; M]; N]) {
    for i in 0..N {
        // println!("Iteration {}", i);
        // matrix.iter().for_each(|line| println!("{:?}", line));
        // println!("");

        // Look for the pivot and swap with the corresponding line
        for k in i..N {
            if matrix[k][i] == 1 {
                matrix.swap(i, k);
                break;
            }
        }

        // Eliminate the other lines
        for k in 0..N {
            if k != i && matrix[k][i] != 0 {
                for j in i..M {
                    matrix[k][j] ^= matrix[i][j];
                }
            }
        }

        // matrix.iter().for_each(|line| println!("{:?}", line));
        // println!("");
    }
}

/// Compute the rank of the given matrix using Gaussian elimination
pub fn rank<const N: usize, const M: usize>(matrix: &[[u8; M]; N]) -> usize {
    // Copy the matrix since the elimination is done in-place
    let mut eliminated = *matrix;
    binary_elimination(&mut eliminated);

    eliminated
        .iter()
        .take_while(|row| row.iter().any(|&coeff| coeff > 0))
        .count()
}

/// Pretty print a matrix represented by an array of lines.
pub fn print_dense<const N: usize, const M: usize>(matrix: &[[u8; M]; N]) {
    matrix.iter().for_each(|line| eprintln!("{:?}", line));
}

/// Count the difference between two slices of the same item type
pub fn count_differences<T: PartialOrd + Copy>(sent: &[T], decoded: &[T]) -> usize {
    sent.iter()
        .zip(decoded)
        .filter(|(&sent_symbol, &decoded_symbol)| sent_symbol != decoded_symbol)
        .count()
}

/// Count the number of combinations of k elements from a n elements set.
pub fn binom(n: usize, k: usize) -> usize {
    if k > n {
        0
    } else {
        (1..=k).fold(1, |acc, val| acc * (n - val + 1) / val)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_difference_count() {
        let count = count_differences(&[1, 4, 5, 6, 5, 3], &[1, 2, 5, 6, 3, 3]);

        assert_eq!(count, 2);
    }

    #[test]
    fn test_rank() {
        #[rustfmt::skip]
        let full_rank = [
            [0, 1, 0, 1],
            [1, 0, 0, 0],
            [1, 1, 0, 0],
        ];
        #[rustfmt::skip]
        let degenerate = [
            [0, 1, 0, 1],
            [1, 0, 0, 1],
            [1, 1, 0, 0],
        ];

        assert_eq!(rank(&full_rank), 3);
        assert_eq!(rank(&degenerate), 2);
    }

    #[test]
    fn test_binom() {
        assert_eq!(binom(1_000, 1), 1_000);
        assert_eq!(binom(1_000, 0), 1);
        assert_eq!(binom(1_000, 4), 41417124750);
    }
}
