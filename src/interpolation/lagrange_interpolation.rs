use crate::interpolation::polynomial::{Polynomial, Evaluatable};


pub fn lagrange_interpolation(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> Polynomial {
    let n = xs.len();

    let mut divided_differences = vec!(ys);

    for i in 1..n {
        let mut di: Vec<f64> = vec!();
        for j in 0..n-i {
            let fx = (divided_differences[i-1][j+1] - divided_differences[i-1][j]) / (xs[j+i] - xs[j]);
            di.push(fx);
        }

        divided_differences.push(di);
    }


    let mut coefficients = vec!();
    for i in 0..n {
        coefficients.push(divided_differences[i][0]);
    }

    Polynomial::new(coefficients)
}


#[cfg(test)]

mod tests {
    use super::*;

    #[test]
    fn test_interpolation_one_point() {
        let xs = vec![1.0];
        let ys = vec![1.0];
        let poly1 = lagrange_interpolation(xs, ys);
        for i in 0..5 {
            assert_eq!(poly1.eval(i as f64), 1.0);
        }
    }
}