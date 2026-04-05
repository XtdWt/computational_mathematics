use pyo3::prelude::*;
use nalgebra::DMatrix;
use crate::interpolation::util::cmp_f64;
use crate::interpolation::polynomial::{Polynomial, Evaluatable, Differentiable, Integrable};


#[pyclass]
pub struct CubicSplineInterpolation {
    pub x_ranges: Vec<(f64, f64)>,
    pub y_functions: Vec<Polynomial>,
}


impl Evaluatable for CubicSplineInterpolation {
    fn eval(&self, x: f64) -> Option<f64> {
        if x < self.x_ranges[0].0 {
            return self.y_functions[0].eval(x)
        }
        for i in 0..self.x_ranges.len() {
            let x_min = self.x_ranges[i].0;
            let x_max = self.x_ranges[i].1;
            if x >= x_min && x <= x_max {
                return self.y_functions[i].eval(x)
            }
        }
        self.y_functions[self.y_functions.len()-1].eval(x)
    }
}


impl Differentiable for CubicSplineInterpolation {
    fn differentiate(&self) -> Self {
        let mut diff_poly = Vec::new();
        for polynomial in self.y_functions.iter() {
            diff_poly.push(polynomial.differentiate())
        }
        CubicSplineInterpolation{
            x_ranges: self.x_ranges.clone(),
            y_functions: diff_poly,
        }
    }
}


impl Integrable for CubicSplineInterpolation {
    fn integrate(&self, x0: f64, y0: f64) -> CubicSplineInterpolation {
        let (i, p) = if x0 < self.x_ranges[0].0 {
            (0, self.y_functions[0].integrate(x0, y0))
        } else if x0 > self.x_ranges[self.x_ranges.len()-1].1 {
            (self.x_ranges.len()-1, self.y_functions[self.y_functions.len()-1].integrate(x0, y0))
        } else {
            let mut i: usize = 0;
            while i < self.x_ranges.len() {
                let x_min = self.x_ranges[i].0;
                let x_max = self.x_ranges[i].1;
                if x0 >= x_min && x0 <= x_max {
                    break
                };
                i += 1;
            }
            (i, self.y_functions[i].integrate(x0, y0))
        };
        println!("{:?}", p);
        let mut integrated_functions = self.y_functions.clone();
        integrated_functions[i] = p.clone();
        println!("integrated functions: {:?}", integrated_functions);

        // integrate from 0 to i
        let mut j = (i-1) as i32;
        let mut x_lower = self.x_ranges[i].0;
        let mut y_lower = p.eval(x_lower).unwrap();
        println!("x_lower: {:?}", x_lower);
        println!("y_lower: {:?}", y_lower);
        while j >= 0 {
            let p_lower = self.y_functions[j as usize].integrate(x_lower, y_lower);
            println!("p_lower: {:?}", p_lower);
            x_lower = self.x_ranges[j as usize].0;
            y_lower = p_lower.eval(x_lower).unwrap();
            integrated_functions[j as usize] = p_lower;
            j -= 1;
        }
        // integrate from i to n
        let mut k = i+1;
        let mut x_upper = self.x_ranges[i].1;
        let mut y_upper = p.eval(x_upper).unwrap();
        println!("x_upper: {:?}", x_upper);
        println!("y_upper: {:?}", y_upper);
        while k < self.x_ranges.len() {
            let p_upper = self.y_functions[k].integrate(x_upper, y_upper);
            println!("p_upper: {:?}", p_upper);
            x_upper = self.x_ranges[k].1;
            y_upper = p_upper.eval(x_upper).unwrap();
            integrated_functions[k] = p_upper;
            k += 1;
        }
        println!("integrated functions: {:?}", integrated_functions);
        CubicSplineInterpolation {
            x_ranges: self.x_ranges.clone(),
            y_functions: integrated_functions,
        }
    }
}


pub fn cubic_spline_interpolation(
    xs: Vec<f64>,
    ys: Vec<f64>
) -> CubicSplineInterpolation {
    let (xs_sorted, ys_sorted): (Vec<f64>, Vec<f64>) = if !xs.is_sorted() {
        let mut joined: Vec<(f64, f64)> = xs.into_iter().zip(ys).collect();
        joined.sort_by(|joined0, joined1| cmp_f64(&joined0.0, &joined1.0));
        joined.into_iter().unzip()
    } else {
        (xs, ys)
    };
    let a_i = ys_sorted;
    let mut h_i = Vec::new();
    for i in 0..(xs_sorted.len()-1) {
        h_i.push(xs_sorted[i+1] - xs_sorted[i]);
    }
    #[allow(non_snake_case)] let mut A = DMatrix::zeros(
        xs_sorted.len(), xs_sorted.len()
    );
    let mut b = DMatrix::zeros(xs_sorted.len(), 1);

    for i in 0..xs_sorted.len() {
        if i == 0 || i == xs_sorted.len()-1 {
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
    for i in 0..(xs_sorted.len()-1) {
        b_i.push(
            (1.0/h_i[i])*(a_i[i+1] - a_i[i]) - ((h_i[i]/3.0)*(2.0*c_i[i] + c_i[i+1]))
        );
        d_i.push(
            (c_i[i+1]-c_i[i])/(3.0*h_i[i])
        );
    }
    let mut x_ranges = Vec::new();
    for i in 1..xs_sorted.len() {
        x_ranges.push((xs_sorted[i-1], xs_sorted[i]));
    }
    let mut y_functions = Vec::new();
    for i in 0..(xs_sorted.len()-1) {
        let weights = vec![a_i[i], b_i[i], c_i[i], d_i[i]];
        let x_i = xs_sorted[i];
        y_functions.push(
            Polynomial{
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
    use crate::interpolation::polynomial::Integrable;
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

    #[test]
    fn test_cubic_spline_interpolation_three_points() {
        let xs = vec![0.0, 1.0, 2.0];
        let ys = vec![5.0, 0.0, 3.0];

        let c = cubic_spline_interpolation(xs, ys);
        assert_eq!(c.eval(0.0).unwrap(), 5.0);
        assert_eq!(c.eval(1.0).unwrap(), 0.0);
        assert_eq!(c.eval(2.0).unwrap(), 3.0);
        assert_eq!(c.eval(0.5).unwrap(), 1.75);
        assert_eq!((c.eval(1.8).unwrap() -  2.016).abs() < 1e-6, true);
    }

    #[test]
    fn test_cubic_spline_interpolation_four_points() {
        let xs = vec![0.0, 1.0, 2.0, 2.5];
        let ys = vec![0.0, 1.0, 8.0, 9.0];

        let c = cubic_spline_interpolation(xs, ys);
        assert_eq!(c.eval(0.0).unwrap(), 0.0);
        assert_eq!(c.eval(1.0).unwrap(), 1.0);
        assert_eq!(c.eval(2.0).unwrap(), 8.0);
        assert_eq!(c.eval(2.5).unwrap(), 9.0);
        assert_eq!((c.eval(0.1).unwrap() - -0.107).abs() < 1e-6, true);
        assert_eq!((c.eval(1.5).unwrap() - 4.60227).abs() < 1e-5, true);
    }

    #[test]
    fn test_cubic_spline_interpolation_sort_points() {
        let xs = vec![0.0, 2.5, 1.0, 2.0];
        let ys = vec![0.0, 9.0, 1.0, 8.0];

        let c = cubic_spline_interpolation(xs, ys);
        assert_eq!(c.eval(0.0).unwrap(), 0.0);
        assert_eq!(c.eval(1.0).unwrap(), 1.0);
        assert_eq!(c.eval(2.0).unwrap(), 8.0);
        assert_eq!(c.eval(2.5).unwrap(), 9.0);
        assert_eq!((c.eval(0.1).unwrap() - -0.107).abs() < 1e-6, true);
        assert_eq!((c.eval(1.5).unwrap() - 4.60227).abs() < 1e-5, true);
    }

    #[test]
    fn test_cubic_spline_interpolation_polynomial_differentiation() {
        let p = Polynomial{
            weights: vec![1.0, 2.0, 3.0],
            x_i: 0.0
        };
        let p_diff = p.differentiate();
        let derivative = |x: f64| 2.0 + 6.0*x;
        for i in -3..3 {
            assert_eq!(p_diff.eval(i as f64).unwrap(), derivative(i as f64));
        }
    }

    #[test]
    fn test_cubic_spline_interpolation_polynomial_integration() {
        let p = Polynomial{
            weights: vec![1.0, 2.0, 3.0],
            x_i: 0.0
        };
        let p_integral = p.integrate(0.0, 1.0);
        let integral = |x: f64| 1.0 + x + x*x + x*x*x;
        for i in -3..3 {
            assert_eq!(p_integral.eval(i as f64).unwrap(), integral(i as f64));
        }
    }

    #[test]
    fn test_cubic_spline_interpolation_spline_differentiation() {
        let p = Polynomial{
            weights: vec![1.0, 2.0, 3.0],
            x_i: 0.0
        };
        let c = CubicSplineInterpolation{
            x_ranges: vec![(-4.0, 4.0)],
            y_functions: vec![p],
        };
        let c_derivative = c.differentiate();
        let derivative = |x: f64| 2.0 + 6.0*x;
        for i in -3..3 {
            assert_eq!(c_derivative.eval(i as f64).unwrap(), derivative(i as f64));
        }
    }

    #[test]
    fn test_cubic_spline_interpolation_spline_integration() {
        let p1 = Polynomial{
            weights: vec![4.0, 2.0],
            x_i: 0.0
        };
        let p2 = Polynomial{
            weights: vec![1.0, 2.0, 3.0],
            x_i: 0.0
        };
        let p3 = Polynomial{
            weights: vec![-3.0, 2.0, 3.0, 4.0],
            x_i: 0.0
        };
        let c = CubicSplineInterpolation{
            x_ranges: vec![(-4.0, -1.0), (-1.0, 1.0), (1.0, 4.0)],
            y_functions: vec![p1, p2, p3],
        };
        let c_integral = c.integrate(0.0, 0.0);
        let weights = vec![
            vec![2.0, 4.0, 1.0],
            vec![0.0, 1.0, 1.0, 1.0],
            vec![3.0, -3.0, 1.0, 1.0, 1.0],
        ];
        for i in 0..3 {
            assert_eq!(c_integral.y_functions[i].weights, weights[i]);
        }
    }

    #[test]
    fn test_cubic_spline_interpolation_extrapolate() {
        let p = Polynomial{
            weights: vec![1.0, 2.0, 3.0],
            x_i: 0.0
        };
        let c = CubicSplineInterpolation{
            x_ranges: vec![(3.0, 4.0)],
            y_functions: vec![p.clone()],
        };

        for i in -3..3 {
            assert_eq!(c.eval(i as f64).unwrap(), p.eval(i as f64).unwrap());
        }
    }
}
