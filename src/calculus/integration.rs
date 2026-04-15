use crate::Function;

pub fn composite_trapezoid_rule (
    f: Function,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> f64 {
    let step = (b-a)/n_buckets as f64;

    let mut total = 0.0;
    for i in 1..n_buckets {
        let x1 = a + i as f64 * step;
        let x2 = a + (i as f64 - 1.0) * step;
        total += step * f(x1) / 2.0;
        total += step * f(x2) / 2.0;
    }
    total
}


pub fn composite_simpsons_rule (
    f: Function,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> f64 {
    let step = (b-a)/n_buckets as f64;

    let mut total = 0.0;
    for i in 1..n_buckets/2 {
        let x1 = a + 2.0 * i as f64 * step;
        let x2 = a + (2.0 * i as f64 - 1.0) * step;
        let x3 = a + (2.0 * i as f64 - 2.0) * step;
        total += step * f(x1) / 3.0;
        total += step * f(x2) / 3.0;
        total += step * f(x3) / 3.0;
    }
    total
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;

    #[test]
    fn integrate_simple_polynomial() {
        let f = Box::new(|x: f64| 2.0*x + 3.0);
        let big_f = |x: f64| x*x + 3.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert!(
            (composite_trapezoid_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );

        assert!(
            (composite_simpsons_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );
    }

    #[test]
    fn integrate_polynomial() {
        let f = Box::new(|x: f64| 5.0*x.powf(4.0) + 3.0*x.powf(2.0) - 14.0*x - 11.0);
        let big_f = |x: f64| x.powf(5.0) + x.powf(3.0) - 7.0*x.powf(2.0) - 11.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert!(
            (composite_trapezoid_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );

        assert!(
            (composite_simpsons_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );
    }

    #[test]
    fn integrate_sine() {
        let f = Box::new(|x: f64| (2.0*x).sin());
        let big_f = |x: f64| -0.5*(2.0*x).cos();
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert!(
            (composite_trapezoid_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );

        assert!(
            (composite_simpsons_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );
    }

    #[test]
    fn integrate_exponential() {
        let f = Box::new(|x: f64| E.powf(-x));
        let big_f = |x: f64| -1.0*E.powf(-x);
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert!(
            (composite_trapezoid_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );

        assert!(
            (composite_simpsons_rule(f.clone(), a, b, 10) - expected_result).abs() > 1e-6
        );
    }
}