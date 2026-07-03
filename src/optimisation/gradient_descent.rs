pub fn gradient_descent<DF: Fn(f64) -> f64>(
    df: &DF,
    x_0: f64,
    step_size: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut x = x_0;
    for _ in 0..n_max {
        let g = df(x);
        if (step_size * g).abs() < eps_tol {
            return x;
        }
        x -= step_size * g;
    }
    return x;
}


#[cfg(test)]
mod tests {
    use crate::assert_almost_eq;
    use std::f64::consts::{E, PI};
    use super::*;

    #[test]
    fn test_gradient_descent_simple() {
        let _f = |x: f64| x.powi(2) + 3.0;
        let df = |x: f64| 2.0*x;

        let result = gradient_descent(&df, 3.0, 0.2, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_gradient_descent_exponential() {
        let _f = |x: f64| -1.0 / E.powf(x.powi(2));
        let df = |x: f64| (2.0 * x) / E.powf(x.powi(2));

        let result = gradient_descent(&df, 2.0, 0.2, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_gradient_descent_sinusoidal() {
        let _f = |x: f64| x.sin();
        let df = |x: f64| x.cos();

        let result = gradient_descent(&df, 1.2, 0.2, 100, 1e-8);

        assert_almost_eq!(result, -PI/2.0);
    }

    #[test]
    fn test_gradient_descent_two_peaks() {
        let _f = |x: f64| x.powi(3) + 3.0*x.powi(2);
        let df = |x: f64| 3.0*x.powi(2) + 6.0*x;

        let result = gradient_descent(&df, 2.0, 0.2, 100, 1e-8);

        assert_almost_eq!(result, 0.0);

        let result = gradient_descent(&df, -1.9, 0.2, 100, 1e-8);

        assert_almost_eq!(result, 0.0);
    }

    #[test]
    fn test_zero_gradient_initial_condition() {
        let _f = |x: f64| x.powi(3) + 3.0*x.powi(2);
        let df = |x: f64| 3.0*x.powi(2) + 6.0*x;

        let result = gradient_descent(&df, -2.0, 0.2, 100, 1e-8);

        assert_almost_eq!(result, -2.0);
    }
}
