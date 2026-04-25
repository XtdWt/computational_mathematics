use crate::interpolation::polynomial::LagrangePolynomial;

pub fn barycentric_lagrange_interpolation(xs: Vec<f64>, ys: Vec<f64>) -> LagrangePolynomial {
    let mut weights = Vec::new();
    let n = xs.len();
    for i in 0..n {
        let mut weight_i = 1.0;
        for k in 0..n {
            if k != i {
                weight_i *= xs[i] - xs[k]
            }
        }
        weights.push(1.0 / weight_i);
    }
    LagrangePolynomial { weights, xs, ys }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::interpolation::polynomial::Evaluatable;

    #[test]
    fn test_interpolation_one_point() {
        let xs = vec![1.0];
        let ys = vec![2.0];
        let poly1 = barycentric_lagrange_interpolation(xs, ys);
        for i in 0..5 {
            assert_eq!(poly1.eval(i as f64).unwrap(), 2.0);
        }
    }

    #[test]
    fn test_interpolation_two_point() {
        let xs = vec![1.0, 2.0];
        let ys = vec![3.0, 5.0];
        let poly1 = barycentric_lagrange_interpolation(xs, ys);
        let expected_fn = |x: f64| 2.0 * x + 1.0;
        for i in 0..5 {
            assert!((poly1.eval(i as f64).unwrap() - expected_fn(i as f64)).abs() < 1e-6);
        }
    }

    #[test]
    fn test_interpolation_three_point() {
        let xs = vec![0.0, 2.0, 3.0];
        let ys = vec![1.0, 2.0, 4.0];
        let poly2 = barycentric_lagrange_interpolation(xs, ys);
        let poly2_y_values = vec![1.0, 1.0, 2.0, 4.0, 7.0];
        for i in 0..5 {
            assert_eq!(
                (poly2.eval(i as f64).unwrap() - poly2_y_values[i]).abs() < 1e-6,
                true
            );
        }
    }

    #[test]
    fn test_interpolation_five_point() {
        let xs = vec![0.0, 1.0, 2.0, 3.0, 4.0];
        let ys = vec![1.0, -1.0, 2.0, 4.0, 3.0];
        let poly2 = barycentric_lagrange_interpolation(xs, ys);
        let poly2_y_values = vec![1.0, -1.0, 2.0, 4.0, 3.0];
        for i in 0..5 {
            assert_eq!(
                (poly2.eval(i as f64).unwrap() - poly2_y_values[i]).abs() < 1e-6,
                true
            );
        }
        assert_eq!((poly2.eval(-1.0).unwrap() - 18.0).abs() < 1e-6, true);
        assert_eq!((poly2.eval(6.0).unwrap() - 4.0).abs() < 1e-6, true);
    }

    #[test]
    fn test_out_of_order_points() {
        let xs = vec![3.0, 2.0, 0.0];
        let ys = vec![4.0, 2.0, 1.0];
        let poly2 = barycentric_lagrange_interpolation(xs, ys);
        let poly2_y_values = vec![1.0, 1.0, 2.0, 4.0, 7.0];
        for i in 0..5 {
            assert_eq!(
                (poly2.eval(i as f64).unwrap() - poly2_y_values[i]).abs() < 1e-6,
                true
            );
        }
    }
}
