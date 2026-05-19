use rayon::prelude::*;

use crate::Function;

pub fn composite_trapezoid_rule (
    f: Function,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> f64 {
    let step = (b-a)/n_buckets as f64;
    return (1..n_buckets+1)
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
    let indices = if n_buckets % 2 == 0 {
        1..n_buckets/2 + 1
    } else {
        1..(n_buckets-1)/2 + 1
    };

    let step = if n_buckets % 2 == 0 {
        (b-a)/n_buckets as f64
    } else {
        (b-a)/(n_buckets - 1) as f64
    };
    return indices
        .into_par_iter()
        .map(|i| {
            let i_64 = i as f64;
            let x0 = a + (2.0 * i_64) * step;
            let x1 = a + (2.0 * i_64 - 1.0) * step;
            let x2 = a + (2.0 * i_64 - 2.0) * step;

            vec![x0, x1, x2]
        }).collect::<Vec<Vec<f64>>>()
        .into_iter()
        .map(|i_list| {(step / 3.0) * (f(i_list[0]) + 4.0 * f(i_list[1]) + f(i_list[2]))})
        .sum();
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;
    use crate::assert_almost_eq;

    #[test]
    fn integrate_simple_polynomial_n_one() {
        let f = Box::new(|x: f64| 2.0*x + 3.0);
        let a = 1.0;
        let b = 4.0;

        assert_almost_eq!(composite_trapezoid_rule(f.clone(), a, b, 1), 24.0);
        assert_almost_eq!(composite_simpsons_rule(f.clone(), a, b, 2), 24.0);
    }

    #[test]
    fn integrate_simple_polynomial() {
        let f = Box::new(|x: f64| 2.0*x + 3.0);
        let big_f = |x: f64| x*x + 3.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 100),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 100),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 101),
            expected_result
        );
    }

    #[test]
    fn integrate_polynomial_n_one() {
        let f = Box::new(|x: f64| 5.0*x.powf(4.0) + 3.0*x.powf(2.0) - 14.0*x - 11.0);
        let a = 1.0;
        let b = 4.0;

        assert_almost_eq!(composite_trapezoid_rule(f.clone(), a, b, 1), 1866.0);
        assert_almost_eq!(composite_simpsons_rule(f.clone(), a, b, 2), 958.125);
    }

    #[test]
    fn integrate_polynomial() {
        let f = Box::new(|x: f64| 5.0*x.powf(4.0) + 3.0*x.powf(2.0) - 14.0*x - 11.0);
        let big_f = |x: f64| x.powf(5.0) + x.powf(3.0) - 7.0*x.powf(2.0) - 11.0*x;
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 10_000),
            expected_result,
            1e-4
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 10_000),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 10_001),
            expected_result
        );
    }

    #[test]
    fn integrate_sine_n_one() {
        let f = Box::new(|x: f64| (2.0*x).sin());
        let a = 1.0;
        let b = 4.0;

        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 1),
            1.5*((8.0_f64).sin() + (2.0_f64).sin())
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 2),
            0.5*((8.0_f64).sin() + 4.0*(5.0_f64).sin() + (2.0_f64).sin())
        );
    }

    #[test]
    fn integrate_sine() {
        let f = Box::new(|x: f64| (2.0*x).sin());
        let big_f = |x: f64| -0.5*(2.0*x).cos();
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 10_000),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 100),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 101),
            expected_result
        );
    }

    #[test]
    fn integrate_exponential_n_one() {
        let f = Box::new(|x: f64| E.powf(-x));
        let a = 1.0;
        let b = 4.0;


        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 1),
            1.5*(E.powi(-4) + E.powi(-1))
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 2),
            0.5*(E.powi(-4) + 4.0*E.powf(-2.5) + E.powi(-1))
        );
    }

    #[test]
    fn integrate_exponential() {
        let f = Box::new(|x: f64| E.powf(-x));
        let big_f = |x: f64| -1.0*E.powf(-x);
        let a = 1.0;
        let b = 4.0;

        let expected_result = big_f(b) - big_f(a);

        assert_almost_eq!(
            composite_trapezoid_rule(f.clone(), a, b, 10_000),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 1_000),
            expected_result
        );

        assert_almost_eq!(
            composite_simpsons_rule(f.clone(), a, b, 1_001),
            expected_result
        );
    }
}
