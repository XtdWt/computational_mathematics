use pyo3::prelude::*;
use crate::interpolation::polynomial::{Evaluatable};


#[pyclass]
pub struct NewtonsDividedDifferencePolynomial {
    divided_differences: Vec<f64>,
    xs: Vec<f64>,
    n: usize,
}


impl Evaluatable for NewtonsDividedDifferencePolynomial {
    fn eval(&self, x: f64) -> f64 {
        let mut total = 0.0;
        for i in 0..self.n {
            let mut c = self.divided_differences[i];
            for j in 0..i {
                c = c * (x - self.xs[j])
            }
            total += c;
        }
        total
    }
}


#[pymethods]
impl NewtonsDividedDifferencePolynomial {
    fn __call__(&self, x:f64) -> f64 {
        self.eval(x)
    }
}


pub fn newtons_divided_difference_interpolation(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> NewtonsDividedDifferencePolynomial {
    let mut divided_differences_table = Vec::new();
    let n = xs.len();
    divided_differences_table.push(ys);
    for i in 1..n {
        let mut di = Vec::new();
        for j in 0..n-i {
            let fx = (divided_differences_table[i-1][j+1] - divided_differences_table[i-1][j]) / (xs[j+i] - xs[j]);
            di.push(fx);
        }

        divided_differences_table.push(di);
    }
    let mut divided_differences = Vec::new();
    for i in 0..n {
        divided_differences.push(divided_differences_table[i][0]);
    }
    NewtonsDividedDifferencePolynomial{
        divided_differences,
        xs,
        n,
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_interpolation_one_point() {
        let xs = vec![1.0];
        let ys = vec![2.0];
        let poly1 = newtons_divided_difference_interpolation(xs, ys);
        for i in 0..5 {
            assert_eq!(poly1.eval(i as f64), 2.0);
        }
    }

    #[test]
    fn test_interpolation_three_point() {
        let xs = vec![0.0, 2.0, 3.0];
        let ys = vec![1.0, 2.0, 4.0];
        let poly2 = newtons_divided_difference_interpolation(xs, ys);
        let poly2_y_values = vec!(1.0, 1.0, 2.0, 4.0, 7.0);
        for i in 0..5 {
            assert_eq!((poly2.eval(i as f64) - poly2_y_values[i]).abs() < 1e-6, true);
        }
    }
}