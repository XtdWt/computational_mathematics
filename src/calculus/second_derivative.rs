use crate::Function;
use crate::calculus::util::DerivativeType;


fn second_order_forward_difference(function: Function, x: f64, h: f64) -> f64 {
    (2.0 * function(x) - 5.0 * function(x + h) + 4.0 * function(x + 2.0 * h)
        - function(x + 3.0 * h))
        / (h.powf(2.0))
}


fn second_order_backward_difference(function: Function, x: f64, h: f64) -> f64 {
    (2.0 * function(x) - 5.0 * function(x - h) + 4.0 * function(x - 2.0 * h)
        - function(x - 3.0 * h))
        / (h.powf(2.0))
}


fn second_order_central_difference(function: Function, x: f64, h: f64) -> f64 {
    (function(x + h) - 2.0 * function(x) + function(x - h)) / (h.powf(2.0))
}


pub fn second_derivative(function: Function, x: f64, h: f64, method: DerivativeType) -> f64 {
    match method {
        DerivativeType::Forward => second_order_forward_difference(function, x, h),
        DerivativeType::Backward => second_order_backward_difference(function, x, h),
        DerivativeType::Central => second_order_central_difference(function, x, h),
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;

    #[test]
    fn second_derivative_simple_polynomial() {
        let f = Box::new(|x: f64| x * x + 3.0);
        let ddf = |_x: f64| 2.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Forward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Backward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Central) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
        }
    }

    #[test]
    fn second_derivative_polynomial() {
        let f = Box::new(
            |x: f64| 5.0 * x.powf(4.0) + 3.0 * x.powf(3.0) - 14.0 * x - 11.0
        );
        let ddf = |x: f64| 60.0 * x.powf(2.0) + 18.0 * x;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Forward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Backward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Central) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
        }
    }

    #[test]
    fn second_derivative_sine() {
        let f = Box::new(|x: f64| (2.0 * x).sin());
        let ddf = |x: f64| -4.0 * ((2.0 * x).sin());

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Forward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Backward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Central) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
        }
    }

    #[test]
    fn second_derivative_exponential() {
        let f = Box::new(|x: f64| E.powf(-x));
        let ddf = |x: f64| E.powf(-x);

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Forward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Backward) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
            assert_eq!(
                ((second_derivative(f.clone(), x, 0.000001, DerivativeType::Central) - ddf(x))
                    / ddf(x))
                .abs()
                    < 1e-2,
                true
            );
        }
    }
}