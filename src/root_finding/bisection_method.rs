use crate::Function;

pub fn bisection_method(f: &Function, a: f64, b: f64, n_max: usize, eps_tol: f64) -> f64 {
    let mut lower = a;
    let mut upper = b;
    let mut f_low = f(lower);
    let mut c = (lower + upper) / 2.0;
    for _ in 0..n_max {
        c = (lower + upper) / 2.0;
        let f_c = f(c);

        if f_c.abs() < eps_tol {
            return c;
        }
        if f_c * f_low < 0.0 {
            upper = c;
        } else {
            lower = c;
            f_low = f_c;
        }
    }
    return c;
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn test_bisection_method_one_iteration() {
        let f: Function = Box::new(|x| x + 2.0);
        assert_eq!(
            bisection_method(&f, -3.0, 0.0, 1, 1e-6),
            -1.5
        );
    }

    #[test]
    fn test_bisection_method_two_iteration() {
        let f: Function = Box::new(|x| x + 2.0);
        assert_eq!(
            bisection_method(&f, -3.0, 0.0, 2, 1e-6),
            -2.25
        );
    }

    #[test]
    fn test_bisection_method_many_iteration() {
        let f: Function = Box::new(|x| x + 2.0);
        assert_almost_eq!(bisection_method(&f, -3.0, 0.0, 20, 1e-6), -2.0, 1e-6);
    }

    #[test]
    fn test_bisection_method_many_iteration_complicated() {
        let f: Function = Box::new(|x| x * x - 1.7);
        assert_almost_eq!(bisection_method(&f, 1.0, 2.0, 20, 1e-6), 1.30384, 1e-6);
    }

    #[test]
    fn test_bisection_method_tolerance_break() {
        let f: Function = Box::new(|x| x * x - 1.7);
        assert_eq!(
            bisection_method(&f, 1.0, 2.0, 20, 1.0),
            bisection_method(&f, 1.0, 2.0, 1, 1e-6)
        );
    }
}
