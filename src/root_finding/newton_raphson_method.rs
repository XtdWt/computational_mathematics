use crate::Function;


pub fn newton_raphson_method(
    function: Function,
    derivative: Function,
    x_0: f64,
    n_max: i64,
) -> f64 {
    let mut x = x_0;

    for _ in 0..n_max {
        let y = function(x);
        let dydx = derivative(x);
        let dx = - y/dydx;
        x = x + dx;
    }
    x
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_newton_raphson_method_one_iteration() {
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        let df = Box::new(|x: f64| 2.0*x + 3.0);
        assert_eq!(
            newton_raphson_method(f, df, 1.0, 1), 1.0 - (1.0/5.0)
        )
    }

    #[test]
    fn test_newton_raphson_method_two_iteration() {
        let epsilon = 1e-6;
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        let df = Box::new(|x: f64| 2.0*x + 3.0);
        assert!(
            newton_raphson_method(f, df, 1.0, 2) - 0.8 - (0.04 / 4.6) < epsilon
        )
    }

    #[test]
    fn test_newton_raphson_method_many_iteration() {
        let epsilon = 1e-6;
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        let df = Box::new(|x: f64| 2.0*x + 3.0);
        assert!(
            newton_raphson_method(f, df, 1.0, 20) - 0.79128784 < epsilon
        )
    }

    #[test]
    fn test_newton_raphson_method_simple() {
        let epsilon = 1e-6;
        let f = Box::new(|x: f64| x + 3.0);
        let df = Box::new(|_x: f64| 1.0);
        assert!(
            newton_raphson_method(f, df, 1.0, 20) + 3.0 < epsilon
        )
    }
}