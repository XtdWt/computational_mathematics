use std::f64::consts::PI;
use crate::interpolation::util::cmp_f64;


pub fn chebyshev_nodes(
    a: f64,
    b: f64,
    n: usize,
) -> Vec<f64> {
    let mut nodes = Vec::new();
    let a1 = (a+b)/2.0;
    let a2 = (b-a)/2.0;
    let n_64 = n as f64;
    for i in 0..n {
        let i_64 = i as f64;
        let a3 = ((PI * (2.0 * i_64 + 1.0))/(2.0 * n_64)).cos();
        nodes.push(a1 + a2 * a3)
    }
    nodes.sort_by(cmp_f64);
    nodes
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
        for i in 0..expected_nodes.len() {
            assert_eq!((nodes[i] - expected_nodes[i]).abs() < 1e-6, true);
        }
    }

    #[test]
    fn chebyshev_nodes_0_1_three_node() {
        let nodes = chebyshev_nodes(0.0, 1.0, 3);
        let expected_nodes = vec![0.0669873, 0.5, 0.9330127];
        for i in 0..expected_nodes.len() {
            assert_eq!((nodes[i] - expected_nodes[i]).abs() < 1e-6, true);
        }
    }

    #[test]
    fn chebyshev_nodes_0_10_five_node() {
        let nodes = chebyshev_nodes(0.0, 10.0, 5);
        let expected_nodes = vec![
            0.24471742,
            2.06107374,
            5.0,
            7.93892626,
            9.75528258,
        ];
        for i in 0..expected_nodes.len() {
            assert_eq!((nodes[i] - expected_nodes[i]).abs() < 1e-6, true);
        }
    }
}