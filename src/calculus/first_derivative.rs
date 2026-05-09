use crate::Function;
use crate::calculus::util::DerivativeType;


fn forward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    return (function(x+h) - function(x)) / h;
}


fn backward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    return (function(x) - function(x-h)) / h;
}


fn central_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    return (function(x + h) - function(x - h)) / (2.0 * h);
}


pub fn first_derivative(
    function: Function,
    x: f64,
    h: f64,
    method: DerivativeType,
) -> f64 {
    return match method {
        DerivativeType::Forward => forward_difference(function, x, h),
        DerivativeType::Backward => backward_difference(function, x, h),
        DerivativeType::Central => central_difference(function, x, h),
    };
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;
    use crate::util::assert_almost_equal;


    #[test]
    fn first_derivative_simple_polynomial() {
        let f = Box::new(|x: f64| x*x + 3.0);
        let df = |x: f64| 2.0*x;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            // assert_eq!(
            //     (
            //         first_derivative(f.clone(), x, 1e-8, DerivativeType::Forward) - df(x)
            //     ).abs() < 1e-6, true
            // );
            assert_almost_equal!(
                first_derivative(f.clone(), x, 1e-8, DerivativeType::Forward),
                df(x)
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Backward) - df(x)
                ).abs() < 1e-6, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Central) - df(x)
                ).abs() < 1e-6, true
            );
        }
    }

    #[test]
    fn first_derivative_polynomial() {
        let f = Box::new(|x: f64| x.powf(3.0) + x.powf(2.0) - 4.0*x - 2.0);
        let df = |x: f64| 3.0*x.powf(2.0) + 2.0*x - 4.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Forward) - df(x)
                ).abs() < 1e-4, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Backward) - df(x)
                ).abs() < 1e-4, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Central) - df(x)
                ).abs() < 1e-4, true
            );
        }
    }

    #[test]
    fn first_derivative_sine() {
        let f = Box::new(|x: f64| (2.0*x).sin());
        let df = |x: f64| 2.0*(2.0*x).cos();

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Forward) - df(x)
                ).abs() < 1e-6, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Backward) - df(x)
                ).abs() < 1e-6, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Central) - df(x)
                ).abs() < 1e-6, true
            );
        }
    }

    #[test]
    fn first_derivative_exponential() {
        let f = Box::new(|x: f64| E.powf(-x));
        let df = |x: f64| -1.0*E.powf(-x);

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Forward) - df(x)
                ).abs() < 1e-6, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Backward) - df(x)
                ).abs() < 1e-6, true
            );
            assert_eq!(
                (
                    first_derivative(f.clone(), x, 1e-8, DerivativeType::Central) - df(x)
                ).abs() < 1e-6, true
            );
        }
    }
}
