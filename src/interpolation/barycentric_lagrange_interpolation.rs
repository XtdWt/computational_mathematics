use pyo3::prelude::*;
use crate::interpolation::polynomial::{Evaluatable};


#[pyclass]
pub struct LagrangePolynomial {
    weights: Vec<f64>,
    xs: Vec<f64>,
    ys: Vec<f64>,
    n: usize,
}


impl Evaluatable for LagrangePolynomial {
    fn eval(&self, x: f64) -> f64 {
        for i in 0..self.n {
            if self.xs[i] == x {
                return self.ys[i];
            }
        }
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for i in 0..self.n {
            numerator += (self.weights[i] * self.ys[i])/(x - self.xs[i]);
            denominator += self.weights[i]/(x - self.xs[i]);
        }
        numerator/denominator
    }
}


#[pymethods]
impl LagrangePolynomial {
    fn __call__(&self, x:f64) -> f64 {
        self.eval(x)
    }
}


pub fn barycentric_lagrange_interpolation(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> LagrangePolynomial {
    let mut weights = Vec::new();
    let n = xs.len();
    for i in 0..n {
        let mut weight_i = 1.0;
        for k in 0..n {
            if k != i {
                weight_i *= xs[i] - xs[k]
            }
        }
        weights.push(1.0/weight_i);
    }
    println!("{:?}", weights);
    LagrangePolynomial{
        weights,
        xs,
        ys,
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
        let poly1 = barycentric_lagrange_interpolation(xs, ys);
        for i in 0..5 {
            assert_eq!(poly1.eval(i as f64), 2.0);
        }
    }

    #[test]
    fn test_interpolation_three_point() {
        let xs = vec![0.0, 2.0, 3.0];
        let ys = vec![1.0, 2.0, 4.0];
        let poly2 = barycentric_lagrange_interpolation(xs, ys);
        let poly2_y_values = vec!(1.0, 1.0, 2.0, 4.0, 7.0);
        for i in 0..5 {
            assert_eq!(poly2.eval(i as f64) - poly2_y_values[i] < 1e-6, true);
        }
    }
}