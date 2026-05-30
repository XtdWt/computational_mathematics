use crate::calculus::composite_integration::composite_trapezoid_rule;


pub fn adaptive_quadrature<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    eps_tol: f64,
) -> f64 {
    let mut ranges = vec![(a, b)];
    let mut total = 0.0;

    while let Some((a, b)) = ranges.pop() {
        let c = 0.5 * (a + b);

        let q  = composite_trapezoid_rule(f, a, b, 1);
        let q1 = composite_trapezoid_rule(f, a, c, 1);
        let q2 = composite_trapezoid_rule(f, c, b, 1);
        if (q - (q1 + q2)).abs() < 3.0 * eps_tol {
            total += q1 + q2;
        } else {
            ranges.push((c, b));
            ranges.push((a, c));
        }
    }
    return total;
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;
    use crate::assert_almost_eq;

    #[test]
    fn integrate_simple_polynomial() {
        let f = |x: f64| 2.0*x + 3.0;
        let big_f = |x: f64| x*x + 3.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(adaptive_quadrature(&f, a, b, 1e-7), expected_result, 1e-5);
    }

    #[test]
    fn integrate_polynomial() {
        let f = |x: f64| 5.0*x.powf(4.0) + 3.0*x.powf(2.0) - 14.0*x - 11.0;
        let big_f = |x: f64| x.powf(5.0) + x.powf(3.0) - 7.0*x.powf(2.0) - 11.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(adaptive_quadrature(&f, a, b, 1e-7), expected_result, 1e-4);
    }

    #[test]
    fn integrate_sine() {
        let f = |x: f64| (2.0*x).sin();
        let big_f = |x: f64| -0.5*(2.0*x).cos();
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(adaptive_quadrature(&f, a, b, 1e-7), expected_result, 1e-5);
    }

    #[test]
    fn integrate_exponential() {
        let f = |x: f64| E.powf(-x);
        let big_f = |x: f64| -1.0*E.powf(-x);
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(adaptive_quadrature(&f, a, b, 1e-7), expected_result, 1e-5);
    }
}
