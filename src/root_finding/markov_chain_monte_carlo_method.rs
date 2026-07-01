use rand::RngExt;
use rand_distr::{Normal, Distribution};

pub fn markov_chain_monte_carlo_method<
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
    C: Fn(f64) -> bool,
>(
    f: &F,
    df: &DF,
    constraint: &C,
    x_0: f64,
    step_size: f64,
    n_max: usize,
    eps_tol: f64,
    sigma: f64,
) -> f64 {
    let mut rng = rand::rng();
    let mut x = x_0;

    for _ in 0..n_max {
        let rand_float: f64 = rng.random();
        let x_next = if rand_float > 0.5 {
            x - step_size * df(x)
        } else {
            let normal = Normal::new(x, sigma).unwrap();
            normal.sample(&mut rng)
        };

        if constraint(x_next) && f(x_next) <= f(x) {
            if (x_next - x).abs() < eps_tol {
                return x;
            }
            x = x_next;
        }
     }
     return x;
}


#[cfg(test)]
mod tests {
    use crate::assert_almost_eq;

use super::*;

    #[test]
    fn test_mcmc_method_quadratic_no_constraint() {
        let f = |x: f64| (x - 3.0).powi(2);
        let df = |x: f64| 2.0 * (x - 3.0);
        let constraint = |_x: f64| true;

        let result = markov_chain_monte_carlo_method(
            &f,
            &df,
            &constraint,
            0.0,
            0.2,
            100,
            1e-8,
            0.05,
        );

        assert_almost_eq!(result, 3.0);
    }

    #[test]
    fn test_mcmc_method_quadratic_constraint() {
        let f = |x: f64| (x - 3.0).powi(2);
        let df = |x: f64| 2.0 * (x - 3.0);
        let constraint = |x: f64| x > 2.5;

        let result = markov_chain_monte_carlo_method(
            &f,
            &df,
            &constraint,
            0.0,
            0.2,
            1000,
            1e-8,
            0.05,
        );

        assert_almost_eq!(result, 2.5);
    }
}
