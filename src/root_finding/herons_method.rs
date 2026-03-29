pub fn herons_method(a: f64, x_0: f64, n_max: i64) -> f64 {
    let mut n: i64 = 0;
    let mut x: f64 = x_0;
    while n < n_max {
        x = 0.5 * (x + a/x);
        n = n + 1;
    }
    x
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_one_iteration() {
        assert_eq!(herons_method(2.0, 1.0, 1), 1.5);
    }

    #[test]
    fn test_two_iteration() {
        assert!(herons_method(2.0, 1.0, 2) - 1.4666666 < 0.0000001);
    }

    #[test]
    fn test_many_iteration() {
        assert!(herons_method(2.0, 1.0, 5) - 2.0_f64.sqrt() < 0.0000001);
    }
}