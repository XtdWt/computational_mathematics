use rayon::prelude::*;

use crate::interpolation::util::cmp_f64;
use std::f64::consts::PI;

pub fn chebyshev_nodes(a: f64, b: f64, n: usize) -> Vec<f64> {
    let a1 = (a + b) / 2.0;
    let a2 = (b - a) / 2.0;
    let n_64 = n as f64;
    let mut nodes: Vec<f64> = (0..n)
        .into_par_iter()
        .map(|i| a1 + a2 * ((PI * (2.0 * i as f64 + 1.0)) / (2.0 * n_64)).cos())
        .collect();
    nodes.sort_by(cmp_f64);
    return nodes;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn chebyshev_nodes_0_1_one_node() {
        let nodes = chebyshev_nodes(0.0, 1.0, 1);
        assert_eq!(nodes, vec![0.5]);
    }

    #[test]
    fn chebyshev_nodes_0_1_two_node() {
        let nodes = chebyshev_nodes(0.0, 1.0, 2);
        let expected_nodes = vec![0.14644661, 0.85355339];
        for (node, expected_node) in nodes.into_iter().zip(expected_nodes) {
            assert_eq!((node - expected_node).abs() < 1e-6, true);
        }
    }

    #[test]
    fn chebyshev_nodes_0_1_three_node() {
        let nodes = chebyshev_nodes(0.0, 1.0, 3);
        let expected_nodes = vec![0.0669873, 0.5, 0.9330127];
        for (node, expected_node) in nodes.into_iter().zip(expected_nodes) {
            assert_eq!((node - expected_node).abs() < 1e-6, true);
        }
    }

    #[test]
    fn chebyshev_nodes_0_10_five_node() {
        let nodes = chebyshev_nodes(0.0, 10.0, 5);
        let expected_nodes = vec![0.24471742, 2.06107374, 5.0, 7.93892626, 9.75528258];
        for (node, expected_node) in nodes.into_iter().zip(expected_nodes) {
            assert_eq!((node - expected_node).abs() < 1e-6, true);
        }
    }
}
