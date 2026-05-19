use crate::Function;
use crate::calculus::util::DerivativeType;


fn second_order_forward_difference(function: &Function, x: f64, h: f64) -> f64 {
    return (function(x + 2.0 * h) - 2.0 * function(x + h) + function(x)) / h.powi(2);
}


fn second_order_backward_difference(function: &Function, x: f64, h: f64) -> f64 {
    return (function(x) - 2.0 * function(x - h) + function(x - 2.0 * h)) / h.powi(2);
}


fn second_order_central_difference(function: &Function, x: f64, h: f64) -> f64 {
    return (function(x + h) - 2.0 * function(x) + function(x - h)) / (h.powi(2));
}


pub fn second_derivative(function: Function, x: f64, h: f64, method: DerivativeType) -> f64 {
    return match method {
        DerivativeType::Forward => second_order_forward_difference(&function, x, h),
        DerivativeType::Backward => second_order_backward_difference(&function, x, h),
        DerivativeType::Central => second_order_central_difference(&function, x, h),
    };
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;
    use crate::assert_almost_eq;

    #[test]
    fn second_derivative_simple_polynomial() {
        let f = Box::new(|x: f64| x * x + 3.0);
        let ddf = |_x: f64| 2.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Forward), ddf(x), 1e-4);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Backward), ddf(x), 1e-4);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Central), ddf(x), 1e-4);
        }
    }

    #[test]
    fn second_derivative_polynomial() {
        let f = Box::new(
            |x: f64| x.powi(3) + 2.0 * x.powi(2) - 14.0 * x - 11.0
        );
        let ddf = |x: f64| 6.0 * x + 4.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Forward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Backward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Central), ddf(x), 1e-3);
        }
    }

    #[test]
    fn second_derivative_sine() {
        let f = Box::new(|x: f64| (2.0 * x).sin());
        let ddf = |x: f64| -4.0 * ((2.0 * x).sin());

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Forward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Backward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Central), ddf(x), 1e-3);
        }
    }

    #[test]
    fn second_derivative_exponential() {
        let f = Box::new(|x: f64| E.powf(-x));
        let ddf = |x: f64| E.powf(-x);

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Forward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Backward), ddf(x), 1e-3);
            assert_almost_eq!(second_derivative(f.clone(), x, 1e-5, DerivativeType::Central), ddf(x), 1e-3);
        }
    }
}
