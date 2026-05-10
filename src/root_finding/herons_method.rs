pub fn herons_method(a: f64, x_0: f64, n_max: usize) -> f64 {
    return (0..n_max).fold(x_0, |total, _| 0.5 * (total + a / total));
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn test_one_iteration() {
        assert_eq!(herons_method(2.0, 1.0, 1), 1.5);
    }

    #[test]
    fn test_two_iteration() {
        assert_almost_eq!(herons_method(2.0, 1.0, 2), 1.4666666, 0.05);
    }

    #[test]
    fn test_many_iteration() {
        assert_almost_eq!(herons_method(2.0, 1.0, 5), 2.0_f64.sqrt());
    }
}
