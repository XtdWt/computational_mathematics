use rayon::prelude::*;

use crate::Function;

pub fn composite_trapezoid_rule (
    f: Function,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> f64 {
    let step = (b-a)/n_buckets as f64;
    return (1..n_buckets)
        .into_par_iter()
        .flat_map(|i| { vec![a + i as f64 * step, a + (i as f64 - 1.0) * step] } )
        .collect::<Vec<f64>>()
        .into_iter()
        .map(|i| {step * f(i) / 2.0})
        .sum();
}


pub fn composite_simpsons_rule (
    f: Function,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> f64 {
    let step = (b-a)/n_buckets as f64;

    return (1..n_buckets/2)
        .into_par_iter()
        .flat_map(|i| { vec![a + 2.0 * i as f64 * step, a + (2.0 * i as f64 - 1.0) * step, a + (2.0 * i as f64 - 2.0) * step] } )
        .collect::<Vec<f64>>()
        .into_iter()
        .map(|i| {step * f(i) / 3.0})
        .sum();
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
