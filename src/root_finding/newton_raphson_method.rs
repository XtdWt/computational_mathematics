use crate::calculus::first_derivative::central_difference;


pub fn quasi_newton_raphson_method<F: Fn(f64) -> f64>(
    f: &F,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut x = x_0;

    for _ in 0..n_max {
        let x_next = x - f(x) / central_difference(&f, x, eps_tol);
        if (x_next - x).abs() < eps_tol {
            return x_next;
        }
        x = x_next;
    }
    return x;
}


pub fn newton_raphson_method<F: Fn(f64) -> f64, DF: Fn(f64) -> f64>(
    f: &F,
    df: &DF,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut x = x_0;

    for _ in 0..n_max {
        let x_next = x - f(x) / df(x);
        if (x_next - x).abs() < eps_tol {
            return x_next;
        }
        x = x_next;
    }
    return x;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn test_newton_raphson_method_one_iteration() {
        let f = |x: f64| x * x + 3.0 * x - 3.0;
        let df = |x: f64| 2.0 * x + 3.0;

        assert_eq!(newton_raphson_method(&f, &df, 1.0, 1, 1e-6), 1.0 - (1.0 / 5.0));
        assert_eq!(quasi_newton_raphson_method(&f, 1.0, 1, 0.1), 1.0 - (1.0 / 5.0));
    }

    #[test]
    fn test_newton_raphson_method_two_iteration() {
        let f = |x: f64| x * x + 3.0 * x - 3.0;
        let df = |x: f64| 2.0 * x + 3.0;
        assert_almost_eq!(newton_raphson_method(&f, &df, 1.0, 2, 1e-6), 0.8 - (0.04 / 4.6));
        assert_almost_eq!(quasi_newton_raphson_method(&f, 1.0, 2, 1e-6), 0.8 - (0.04 / 4.6));
    }

    #[test]
    fn test_newton_raphson_method_many_iteration() {
        let f = |x: f64| x * x + 3.0 * x - 3.0;
        let df = |x: f64| 2.0 * x + 3.0;
        assert_almost_eq!(newton_raphson_method(&f, &df, 1.0, 20, 1e-6), 0.79128784);
        assert_almost_eq!(quasi_newton_raphson_method(&f, 1.0, 20, 1e-6), 0.79128784);
    }

    #[test]
    fn test_newton_raphson_method_simple() {
        let f = |x: f64| x + 3.0;
        let df = |_x: f64| 1.0;
        assert_almost_eq!(newton_raphson_method(&f, &df, 1.0, 20, 1e-6), -3.0);
        assert_almost_eq!(quasi_newton_raphson_method(&f, 1.0, 20, 1e-6), -3.0);
    }

    #[test]
    fn test_newton_raphson_method_tolerance_break() {
        let f = |x: f64| x + 3.0;
        let df = |_x: f64| 1.0;
        assert_almost_eq!(
            newton_raphson_method(&f, &df, 1.0, 20, 1.0),
            newton_raphson_method(&f, &df, 1.0, 1, 1e-6)
        );
        assert_almost_eq!(
            quasi_newton_raphson_method(&f, 1.0, 20, 1.0),
            quasi_newton_raphson_method(&f, 1.0, 1, 1e-6)
        );
    }
}
