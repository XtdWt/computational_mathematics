use crate::calculus::util::DerivativeType;


fn forward_difference<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    h: f64,
) -> f64 {
    return (f(x+h) - f(x)) / h;
}


fn backward_difference<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    h: f64,
) -> f64 {
    return (f(x) - f(x-h)) / h;
}


fn central_difference<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    h: f64,
) -> f64 {
    return (f(x + h) - f(x - h)) / (2.0 * h);
}


pub fn first_derivative<F: Fn(f64) -> f64>(
    f: &F,
    x: f64,
    h: f64,
    method: DerivativeType,
) -> f64 {
    return match method {
        DerivativeType::Forward => forward_difference(f, x, h),
        DerivativeType::Backward => backward_difference(f, x, h),
        DerivativeType::Central => central_difference(f, x, h),
    };
}


#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::E;
    use crate::assert_almost_eq;

    #[test]
    fn first_derivative_simple_polynomial_h_one() {
        let f = |x: f64| x*x + 3.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Forward),
                2.0*x+1.0
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Backward),
                2.0*x-1.0
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Central),
                2.0*x
            );
        }
    }

    #[test]
    fn first_derivative_simple_polynomial() {
        let f = |x: f64| x*x + 3.0;
        let df = |x: f64| 2.0*x;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Forward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Backward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Central),
                df(x),
                1e-6
            );
        }
    }

    #[test]
    fn first_derivative_polynomial_h_one() {
        let f = |x: f64| x.powf(3.0) + x.powf(2.0) - 4.0*x - 2.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Forward),
                3.0*x.powi(2) + 5.0*x - 2.0
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Backward),
                3.0*x.powi(2) - x - 4.0
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Central),
                3.0*x.powi(2) + 2.0*x - 3.0
            );
        }
    }

    #[test]
    fn first_derivative_polynomial() {
        let f = |x: f64| x.powf(3.0) + x.powf(2.0) - 4.0*x - 2.0;
        let df = |x: f64| 3.0*x.powf(2.0) + 2.0*x - 4.0;

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Forward),
                df(x),
                1e-4
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Backward),
                df(x),
                1e-4
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Central),
                df(x),
                1e-4
            );
        }
    }

    #[test]
    fn first_derivative_sine_h_one() {
        let f = |x: f64| (2.0*x).sin();

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Forward),
                (2.0*x + 2.0).sin() - (2.0*x).sin()
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Backward),
                (2.0*x).sin() - (2.0*x - 2.0).sin()
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Central),
                ((2.0*x + 2.0).sin() - (2.0*x - 2.0).sin())/2.0
            );
        }
    }

    #[test]
    fn first_derivative_sine() {
        let f = |x: f64| (2.0*x).sin();
        let df = |x: f64| 2.0*(2.0*x).cos();

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Forward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Backward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Central),
                df(x),
                1e-6
            );
        }
    }

    #[test]
    fn first_derivative_exponential_h_one() {
        let f = |x: f64| E.powf(-x);

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Forward),
                E.powf(-x)*(E.powi(-1) - 1.0)
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Backward),
                E.powf(-x)*(1.0 - E)
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1.0, DerivativeType::Central),
                E.powf(-x)*(E.powi(-1) - E.powi(1)) / 2.0
            );
        }
    }

    #[test]
    fn first_derivative_exponential() {
        let f = |x: f64| E.powf(-x);
        let df = |x: f64| -1.0*E.powf(-x);

        let xs = vec![-3.3, -2.8, -1.2, -0.7, 0.4, 1.5, 2.9, 3.4];
        for x in xs {
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Forward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Backward),
                df(x),
                1e-6
            );
            assert_almost_eq!(
                first_derivative(&f, x, 1e-8, DerivativeType::Central),
                df(x),
                1e-6
            );
        }
    }
}
