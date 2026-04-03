use pyo3::prelude::*;
use nalgebra::DMatrix;
use crate::interpolation::util::cmp_f64;
use crate::interpolation::util::Evaluatable;


struct CubicPolynomial {
    weights: [f64; 4],
    x_i: f64,
}


impl Evaluatable for CubicPolynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        let mut total = 0.0;
        let x = x - self.x_i;
        let mut x_power = 1.0;
        for w in &self.weights {
            total += w * x_power;
            x_power *= x;
        }
        Some(total)
    }
}


#[pyclass]
pub struct CubicSplineInterpolation {
    x_ranges: Vec<(f64, f64)>,
    y_functions: Vec<CubicPolynomial>,
}


impl Evaluatable for CubicSplineInterpolation {
    fn eval(&self, x: f64) -> Option<f64> {
        for i in 0..self.x_ranges.len() {
            let x_min = self.x_ranges[i].0;
            let x_max = self.x_ranges[i].1;
            if x >= x_min && x <= x_max {
                return self.y_functions[i].eval(x)
            }
        }
        None
    }
}


pub fn cubic_spline_interpolation(
    mut xs: Vec<f64>,
    ys: Vec<f64>
) -> CubicSplineInterpolation {
    let a_i = ys;
    let mut h_i = Vec::new();
    for i in 0..(xs.len()-1) {
        h_i.push(xs[i+1] - xs[i]);
    }
    let mut A = DMatrix::zeros(xs.len(), xs.len());
    let mut b = DMatrix::zeros(xs.len(), 1);

    for i in 0..xs.len() {
        if i == 0 || i == xs.len()-1 {
            A[(i, i)] = 1.0;
            b[(i, 0)] = 0.0;

        } else {
            A[(i, i-1)] = h_i[i-1];
            A[(i, i)] = (h_i[i-1] + h_i[i])*2.0;
            A[(i, i+1)] = h_i[i];
            b[(i, 0)] = (3.0/h_i[i]) * (a_i[i+1] - a_i[i]) - (3.0/h_i[i-1]) * (a_i[i] - a_i[i-1]);
        }
    }
    let c_i: Vec<f64> = (A.try_inverse().unwrap() * b).column(0).iter().cloned().collect();
    let mut b_i = Vec::new();
    let mut d_i = Vec::new();
    for i in 0..(xs.len()-1) {
        b_i.push(
            (1.0/h_i[i])*(a_i[i+1] - a_i[i]) - ((h_i[i]/3.0)*(2.0*c_i[i] + c_i[i+1]))
        );
        d_i.push(
            (c_i[i+1]-c_i[i])/(3.0*h_i[i])
        );
    }
    let mut x_ranges = Vec::new();
    xs.sort_by(cmp_f64);
    for i in 1..xs.len() {
        x_ranges.push((xs[i-1], xs[i]));
    }
    let mut y_functions = Vec::new();
    for i in 0..(xs.len()-1) {
        let weights = [a_i[i], b_i[i], c_i[i], d_i[i]];
        let x_i = xs[i];
        y_functions.push(
            CubicPolynomial{
                weights,
                x_i
            }
        )
    }
    CubicSplineInterpolation{
        x_ranges,
        y_functions,
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_cubic_spline_interpolation_two_points() {
        let xs = vec![0.0, 1.0];
        let ys = vec![0.0, 1.0];

        let c = cubic_spline_interpolation(xs, ys);
        assert_eq!(c.eval(0.0).unwrap(), 0.0);
        assert_eq!(c.eval(1.0).unwrap(), 1.0);
        assert_eq!(c.eval(0.5).unwrap(), 0.5);
    }
}