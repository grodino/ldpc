//! Decoding of error correcting codes

use std::collections::HashMap;
use std::iter;

use crate::channel::{BinaryErasureChannel, BinaryMemorylessChannel, Channel};
use crate::coding::{Code, Ldpc};

pub trait Decoder<'a, CH: Channel, CO: Code<N>, const N: usize> {
    // type Output: Copy + PartialOrd + Sync + Default;

    fn decode(&self, received: &[CH::Output; N]) -> [CO::Output; N];
    fn new(channel: &CH, code: &'a CO) -> Self;
}

#[derive(Debug)]
struct Message<T> {
    bit_to_check: T,
    check_to_bit: T,
}

pub struct BECDecoder<'a, const N: usize, const M: usize> {
    pub erasure_proba: f64,
    pub code: &'a Ldpc<N, M>,
}

struct BECDecoderState<const N: usize> {
    initial_belief: [u8; N],
    belief: [u8; N],
    messages: HashMap<(usize, usize), Message<u8>>,
}

impl<const N: usize, const M: usize> BECDecoder<'_, N, M> {
    fn init_state(&self, received: &[u8; N]) -> BECDecoderState<N> {
        // Store the two messages associated to each vertex (one bit to node
        // message and one node to bit) in a map. NOTE: The hash map could be
        // better choosen since here we have a perfect hash already (```bit_node
        // * n_bit_nodes + check_node```) TODO: create a MessageStore structure
        // deriving std::ops::IndexMut
        let mut messages = HashMap::with_capacity(self.code.total_degree as usize);
        let initial_belief = *received;

        for (bit_node, check_neighbors) in self.code.parity_col.iter().enumerate() {
            for &check_node in check_neighbors.iter() {
                messages.insert(
                    (bit_node, check_node),
                    Message {
                        bit_to_check: initial_belief[bit_node], // Initial bit to check messages are the received symbols
                        check_to_bit: 2, // Initial check to bit messages are erasures
                    },
                );
            }
            // Add the pseudo check-node that always sends the received bits as
            // message
            messages.insert(
                (bit_node, M),
                Message {
                    bit_to_check: 0,
                    check_to_bit: initial_belief[bit_node],
                },
            );
        }

        BECDecoderState {
            messages,
            initial_belief,
            belief: initial_belief,
        }
    }

    /// Bit to check : iterate over the bit nodes to compute the message they
    /// will send to each check node
    fn bit_to_check(&self, state: &mut BECDecoderState<N>) {
        self.code
            .parity_col
            .iter()
            .enumerate()
            .for_each(|(bit_node, neighbor_checknodes)| {
                // println!("\nbit node {bit_node}");

                neighbor_checknodes.iter().for_each(|&check_node_dest| {
                    // Sum the incoming messages (incoming are check nodes AND
                    // message's bits !)
                    let out = match neighbor_checknodes
                        .iter()
                        .chain(iter::once(&M)) // Add the pseudo checknode representing the message
                        // Do not take the current edge (bit_node, check_node_dest) into account
                        .filter(|&&check_node_src| check_node_src != check_node_dest)
                        .fold((true, 2), |(all_erasures, accum), &check_node_src| {
                            let check2bit =
                                state.messages[&(bit_node, check_node_src)].check_to_bit;
                            // let check2bit = check_to_bit[check_node_src][bit_node];
                            let message = match check2bit {
                                2 => accum,
                                m => m,
                            };

                            // For debug only: if incoming messages are not
                            // all erasures then they must agree
                            if !all_erasures && check2bit != 2 && message != check2bit {
                                panic!("AHHHHH{check_node_src:?} {accum}");
                            }

                            (all_erasures & (check2bit == 2), message)
                        }) {
                        (false, a) => a,
                        (true, _) => 2,
                    };

                    // println!("{bit_node}--[{out}]-->{}", check_node_dest);
                    state
                        .messages
                        .get_mut(&(bit_node, check_node_dest))
                        .unwrap()
                        .bit_to_check = out
                });
            });
    }

    /// Check to bit : iterate over the check nodes to compute the message they
    /// will send to each bit node
    fn check_to_bit(&self, state: &mut BECDecoderState<N>) {
        self.code
            .parity_row
            .iter()
            .enumerate()
            .for_each(|(check_node, neighbor_bitnodes)| {
                // println!("\ncheck node {check_node}");

                neighbor_bitnodes.iter().for_each(|&bit_node_dest| {
                    // Sum the messages coming from linked bit nodes
                    let mut out = Some(0);

                    for &bit_node_src in neighbor_bitnodes
                        .iter()
                        .filter(|&&src_bit| src_bit != bit_node_dest)
                    {
                        let bit2check = state.messages[&(bit_node_src, check_node)].bit_to_check;

                        if bit2check == 2 {
                            // If erasure
                            out = None;
                            break;
                        } else {
                            out = Some((out.unwrap() + bit2check) % 2);
                        }
                    }

                    // println!("{bit_node_dest}<--[{out}]--{check_node}");
                    // If out is None, it means that there was an erasure
                    state
                        .messages
                        .get_mut(&(bit_node_dest, check_node))
                        .unwrap()
                        .check_to_bit = out.unwrap_or(2);
                });
            });
    }

    /// Compute the beliefs at each bit node
    fn belief(&self, state: &mut BECDecoderState<N>) {
        for (bit_node, neighbor_checknodes) in self.code.parity_col.iter().enumerate() {
            if state.initial_belief[bit_node] != 2 {
                state.belief[bit_node] = state.initial_belief[bit_node];
            } else {
                state.belief[bit_node] = match neighbor_checknodes.iter().fold(
                    (true, 2),
                    |(all_erasures, accum), &check_node_src| {
                        let check2bit = state.messages[&(bit_node, check_node_src)].check_to_bit;
                        let message = match check2bit {
                            2 => accum,
                            m => m,
                        };

                        // For debug only: if incoming messages are not
                        // all erasures then they must agree
                        if !all_erasures && check2bit != 2 && message != check2bit {
                            panic!("AHHHHH{check_node_src:?} {accum}");
                        }

                        (all_erasures & (check2bit == 2), message)
                    },
                ) {
                    (false, a) => a,
                    (true, _) => 2,
                };
            }
        }
    }
}

impl<'a, const M: usize, const N: usize> Decoder<'a, BinaryErasureChannel, Ldpc<N, M>, N>
    for BECDecoder<'a, N, M>
{
    // type Output = u8;

    fn decode(&self, received: &[u8; N]) -> [u8; N] {
        let state = &mut self.init_state(received);
        let mut i = 0;

        while state.belief.contains(&2) && i < 10 {
            // println!("\n Iteration {i}");
            self.bit_to_check(state);

            self.belief(state);
            // println!("belief {belief:?}");

            self.check_to_bit(state);

            i += 1;
        }

        state.belief
    }

    fn new(channel: &BinaryErasureChannel, code: &'a Ldpc<N, M>) -> Self {
        Self {
            erasure_proba: channel.erasure_proba,
            code,
        }
    }
}

// TODO : create new() method
pub struct BMCDecoder<'a, const N: usize, const M: usize> {
    pub noise_variance: f64,
    pub code: &'a Ldpc<N, M>,
}

struct BMCDecoderState<const N: usize> {
    initial_belief: [f64; N],
    belief: [f64; N],
    messages: HashMap<(usize, usize), Message<f64>>,
}

impl<const N: usize, const M: usize> BMCDecoder<'_, N, M> {
    fn init_state(&self, received: &[f64; N]) -> BMCDecoderState<N> {
        // Store the two messages associated to each vertex (one bit to node
        // message and one node to bit) in a map. NOTE: The hash map could be
        // better choosen since here we have a perfect hash already (```bit_node
        // * n_bit_nodes + check_node```) TODO: create a MessageStore structure
        // deriving std::ops::IndexMut
        let mut messages = HashMap::with_capacity(self.code.total_degree as usize);

        // Compute the initial belief. TODO : replace this by the channel belief
        // function. NOTE: array.map is said not very efficient by rust doc on
        // large arrays.
        let initial_belief = received.map(|amplitude| 2.0 * amplitude / self.noise_variance);

        for (bit_node, check_neighbors) in self.code.parity_col.iter().enumerate() {
            for &check_node in check_neighbors.iter() {
                messages.insert(
                    (bit_node, check_node),
                    Message {
                        // Initial bit to check messages are the received symbols
                        bit_to_check: initial_belief[bit_node],
                        // Initial check to bit messages are erasures
                        check_to_bit: 0.0,
                    },
                );
            }
        }

        BMCDecoderState {
            initial_belief,
            belief: initial_belief,
            messages,
        }
    }

    /// Bit to check : iterate over the bit nodes to compute the message they
    /// will send to each check node
    fn bit_to_check(&self, state: &mut BMCDecoderState<N>) {
        self.code
            .parity_col
            .iter()
            .enumerate()
            .for_each(|(bit_node, neighbor_checknodes)| {
                // println!("\nbit node {bit_node}");

                for &check_node_dest in neighbor_checknodes {
                    let check2bit = state.messages[&(bit_node, check_node_dest)].check_to_bit;

                    // println!("{bit_node}--[{out}]-->{}", check_node_dest);
                    state
                        .messages
                        .get_mut(&(bit_node, check_node_dest))
                        .unwrap()
                        .bit_to_check = state.belief[bit_node] - check2bit;
                }
            });
    }

    /// Check to bit : iterate over the check nodes to compute the message they
    /// will send to each bit node
    fn check_to_bit(&self, state: &mut BMCDecoderState<N>) {
        // Check to bit : enumerate over check node to compute the message they will
        // send to bit nodes
        self.code
            .parity_row
            .iter()
            .enumerate()
            .for_each(|(check_node_src, neighbor_bitnodes)| {
                // println!("\ncheck node {check_node}");

                for &bit_node_dest in neighbor_bitnodes {
                    // Compute the outgoing message from the incoming messages from
                    // linked bit nodes
                    let result = 2.0
                        * neighbor_bitnodes
                            .iter()
                            .filter(|&&src_bit| src_bit != bit_node_dest)
                            .map(|&bit_node| {
                                (state.messages[&(bit_node, check_node_src)].bit_to_check / 2.0f64)
                                    .tanh()
                            })
                            .product::<f64>()
                            .atanh();

                    // println!("{bit_node_dest}<--[{out}]--{check_node}");
                    state
                        .messages
                        .get_mut(&(bit_node_dest, check_node_src))
                        .unwrap()
                        .check_to_bit = result;
                }
            });
    }

    /// Compute the beliefs at each bit node
    fn belief(&self, state: &mut BMCDecoderState<N>) {
        for (bit_node, neighbor_checknodes) in self.code.parity_col.iter().enumerate() {
            state.belief[bit_node] = state.initial_belief[bit_node]
                + neighbor_checknodes
                    .iter()
                    .map(|&check_node| state.messages[&(bit_node, check_node)].check_to_bit)
                    .sum::<f64>()
        }
    }
}

impl<'a, const M: usize, const N: usize> Decoder<'a, BinaryMemorylessChannel, Ldpc<N, M>, N>
    for BMCDecoder<'a, N, M>
{
    // type Output = u8;

    fn decode(&self, received: &[f64; N]) -> [u8; N] {
        let state = &mut self.init_state(received);
        let mut i = 0;

        while !check_syndrome(&state.belief, &self.code.parity_check) && i < 100 {
            // println!("\n Iteration {i}");
            self.check_to_bit(state);

            self.belief(state);
            // println!("belief {:?}", state.belief);

            self.bit_to_check(state);

            i += 1;
            // println!("{}", check_syndrome(&state.belief, &self.code.parity_check));
        }

        let mut decoded = [0u8; N];
        for (i, &ratio) in state.belief.iter().enumerate() {
            decoded[i] = if ratio > 0.0 { 0 } else { 1 };
        }

        decoded
    }

    fn new(channel: &BinaryMemorylessChannel, code: &'a Ldpc<N, M>) -> BMCDecoder<'a, N, M> {
        Self {
            noise_variance: channel.noise_variance,
            code,
        }
    }
}

fn check_syndrome<const N: usize, const M: usize>(
    log_likelyhood_ratio: &[f64; N],
    parity_check: &[[u8; N]; M],
) -> bool {
    // Build the estimated codeword
    let mut estimated = [0u8; N];
    for (i, &ratio) in log_likelyhood_ratio.iter().enumerate() {
        estimated[i] = if ratio >= 0.0 { 0 } else { 1 };
    }

    // Compute and check the syndrome
    for check_equation in parity_check {
        match check_equation
            .iter()
            .zip(estimated.iter())
            .map(|(&estimated_bit, &check_bit)| check_bit * estimated_bit)
            .sum::<u8>()
            % 2
        {
            0 => {}
            _ => return false,
        }
    }

    true
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn bec_decoding() {
        // Length of the code. Also number of bit nodes
        const N: usize = 6;
        // M = N - K the remainder of the code. Also number of check nodes
        const M: usize = 4;
        let parity_check: [[u8; N]; M] = [
            [1, 1, 0, 1, 0, 0],
            [0, 1, 1, 0, 1, 0],
            [1, 0, 0, 0, 1, 1],
            [0, 0, 1, 1, 0, 1],
        ];

        let ldpc = Ldpc::from_dense_parity_check(parity_check);
        let bec_channel = BinaryErasureChannel { erasure_proba: 0.1 };
        let decoder = BECDecoder {
            erasure_proba: bec_channel.erasure_proba,
            code: &ldpc,
        };

        let sent = [0, 0, 1, 0, 1, 1];
        let received: [u8; N] = [0, 0, 1, 2, 2, 2];

        assert_eq!(decoder.decode(&received), sent);
    }

    #[test]
    fn bmc_decoding() {
        // Length of the code. Also number of bit nodes
        const N: usize = 6;
        // M = N - K the remainder of the code. Also number of check nodes
        const M: usize = 4;
        let parity_check: [[u8; N]; M] = [
            [1, 1, 0, 1, 0, 0],
            [0, 1, 1, 0, 1, 0],
            [1, 0, 0, 0, 1, 1],
            [0, 0, 1, 1, 0, 1],
        ];

        let ldpc = Ldpc::from_dense_parity_check(parity_check);
        let bmc_channel = BinaryMemorylessChannel {
            noise_variance: 0.1,
        };
        let decoder = BMCDecoder {
            noise_variance: bmc_channel.noise_variance,
            code: &ldpc,
        };

        let sent = [0, 0, 1, 0, 1, 1];
        let received = [-0.1, 0.5, -0.8, 1.0, -0.7, 0.5];

        assert_eq!(decoder.decode(&received), sent);
    }
}
