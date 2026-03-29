use crate::F;

pub fn bisection_method(
    function: F,
    a: f64,
    b: f64,
    n_max: i64
) -> Option<f64> {
    if function(a) * function(b) > 0.0 {
        return None
    }
    let mut n = 0;
    let mut lower = a;
    let mut upper = b;
    let mut c = (lower + upper) / 2.0;
    while n < n_max {
        c = (lower + upper) / 2.0;
        n = n + 1;
        if function(c) * function(lower) < 0.0 {
            upper = c;
        } else {
            lower = c;
        }
    }
    Some(c)
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bisection_method_one_iteration() {
        assert_eq!(bisection_method(Box::new(|x| x + 2.0), -3.0, 0.0, 1), Some(-1.5));
    }

    #[test]
    fn test_bisection_method_two_iteration() {
        assert_eq!(bisection_method(Box::new(|x| x + 2.0), -3.0, 0.0, 2), Some(-2.25));
    }

    #[test]
    fn test_bisection_method_many_iteration() {
        let epsilon = 1e-6;
        assert!(
            bisection_method(Box::new(|x| x + 2.0), -3.0, 0.0, 20)
                .unwrap() + 2.0 < epsilon
        );
    }

    #[test]
    fn test_bisection_method_many_iteration_complicated() {
        let epsilon = 1e-6;
        assert!(
            bisection_method(Box::new(|x| x * x - 1.7), 1.0, 2.0, 20)
                .unwrap() - 1.30384 < epsilon
        );
    }

    #[test]
    fn test_bisection_method_fails() {
        assert_eq!(bisection_method(Box::new(|x| x + 2.0), -3.0, -4.0, 1), None);
    }
}