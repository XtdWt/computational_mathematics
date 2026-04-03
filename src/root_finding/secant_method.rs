use crate::Function;


pub fn secant_method(
    f: Function,
    x_0: f64,
    x_1: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut x0 = x_0;
    let mut x1 = x_1;
    for _ in 0..n_max {
        let x_next = x1 - f(x1) * (x1 - x0)/(f(x1) - f(x0));
        x0 = x1;
        if (x_next - x1).abs() < eps_tol {
            return x_next;
        }
        x1 = x_next;
    }
    x1
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_secant_method_one_iteration() {
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        assert_eq!(
            secant_method(f, 1.0, 0.9, 1, 1e-6), 0.9 - 0.51*(0.1/0.49)
        )
    }

    #[test]
    fn test_secant_method_two_iteration() {
        let x_2 = 0.9 - 0.51*(0.1/0.49);
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        let f_x_2 = f(x_2);
        let f_x_1 = f(0.9);

        assert_eq!(
            secant_method(
                f, 1.0, 0.9, 2, 1e-6
            ), x_2 - f_x_2*((x_2 - 0.9)/ (f_x_2 - f_x_1))
        )
    }

    #[test]
    fn test_secant_method_many_iteration() {
        let epsilon = 1e-6;
        let f = Box::new(|x: f64| x*x + 3.0*x - 3.0);
        assert!(
            secant_method(f, 1.0, 0.5, 5, 1e-6) - 0.79128784 < epsilon
        )
    }
}