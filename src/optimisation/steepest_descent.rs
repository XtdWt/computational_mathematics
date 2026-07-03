use crate::optimisation::golden_section_search::golden_section_search;

pub fn steepest_descent<
    F: Fn(f64) -> f64,
    DF: Fn(f64) -> f64,
>(
    f: &F,
    df: &DF,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut x = x_0;
    for _ in 0..n_max {
        let g = df(x);
        let line_search_fn = |alpha| f(x - alpha * g);
        let step_size = golden_section_search(&line_search_fn, 0.0, 2.0, 100, 1e-6);
        if (step_size * g).abs() < eps_tol {
            return x;
        }
        x = x - step_size * g;
    }
    return x;
}


#[cfg(test)]
mod tests {
    use crate::assert_almost_eq;
    use std::f64::consts::{E, PI};
    use super::*;

    #[test]
    fn test_steepest_descent_simple() {
        let f = |x: f64| x.powi(2) + 3.0;
        let df = |x: f64| 2.0*x;

        let result = steepest_descent(&f, &df, 3.0, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_steepest_descent_exponential() {
        let f = |x: f64| -1.0 / E.powf(x.powi(2));
        let df = |x: f64| (2.0 * x) / E.powf(x.powi(2));

        let result = steepest_descent(&f, &df, 2.0, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_steepest_descent_sinusoidal() {
        let f = |x: f64| x.sin();
        let df = |x: f64| x.cos();

        let result = steepest_descent(&f, &df, 1.2, 100, 1e-8);

        assert_almost_eq!(result, -PI/2.0);
    }

    #[test]
    fn test_steepest_descent_two_peaks() {
        let f = |x: f64| x.powi(3) + 3.0*x.powi(2);
        let df = |x: f64| 3.0*x.powi(2) + 6.0*x;

        let result = steepest_descent(&f, &df, 2.0, 100, 1e-8);

        assert_almost_eq!(result, 0.0);

        let result = steepest_descent(&f, &df, -1.9, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_zero_gradient_initial_condition() {
        let f = |x: f64| x.powi(3) + 3.0*x.powi(2);
        let df = |x: f64| 3.0*x.powi(2) + 6.0*x;

        let result = steepest_descent(&f, &df, -2.0, 100, 1e-8);

        assert_almost_eq!(result, -2.0);
    }
}
