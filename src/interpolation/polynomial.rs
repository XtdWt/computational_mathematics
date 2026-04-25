use pyo3::pyclass;
use rayon::prelude::*;

pub trait Evaluatable {
    fn eval(&self, x: f64) -> Option<f64>;
}

pub trait Differentiable {
    fn differentiate(&self) -> Self;
}

pub trait Integrable {
    fn integrate(&self, x0: f64, y0: f64) -> Self;
}

#[pyclass]
#[derive(Clone, Debug)]
pub struct Polynomial {
    pub weights: Vec<f64>,
    pub x_i: f64,
}

impl Evaluatable for Polynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        let f_x: f64 = self
            .weights
            .par_iter()
            .enumerate()
            .map(|(i, &w_i)| w_i * (x - self.x_i).powi(i as i32))
            .sum();
        return Some(f_x);
    }
}

impl Differentiable for Polynomial {
    fn differentiate(&self) -> Self {
        let mut diff_weights = Vec::new();
        for i in 1..self.weights.len() {
            diff_weights.push(self.weights[i] * i as f64);
        }

        return Polynomial {
            weights: diff_weights,
            x_i: self.x_i,
        };
    }
}

impl Integrable for Polynomial {
    fn integrate(&self, x0: f64, y0: f64) -> Self {
        let mut int_weights = Vec::new();
        int_weights.push(0.0);
        for i in 0..self.weights.len() {
            int_weights.push(self.weights[i] / (i + 1) as f64);
        }

        let mut definite_integral = Polynomial {
            weights: int_weights,
            x_i: self.x_i,
        };
        let w_0 = y0 - definite_integral.eval(x0).unwrap();
        definite_integral.weights[0] = w_0;
        return definite_integral;
    }
}

#[pyclass]
pub struct NewtonsDividedDifferencePolynomial {
    pub divided_differences: Vec<f64>,
    pub xs: Vec<f64>,
}

impl Evaluatable for NewtonsDividedDifferencePolynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        let mut total = 0.0;
        for i in 0..self.xs.len() {
            let mut c = self.divided_differences[i];
            for j in 0..i {
                c = c * (x - self.xs[j])
            }
            total += c;
        }
        return Some(total);
    }
}

#[pyclass]
pub struct LagrangePolynomial {
    pub weights: Vec<f64>,
    pub xs: Vec<f64>,
    pub ys: Vec<f64>,
}

impl Evaluatable for LagrangePolynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        for i in 0..self.xs.len() {
            if self.xs[i] == x {
                return Some(self.ys[i]);
            }
        }

        let (numerator, denominator) = self
            .xs
            .par_iter()
            .zip(&self.ys)
            .zip(&self.weights)
            .map(|((&xi, &yi), &wi)| {
                let diff = x - xi;
                (wi * yi / diff, wi / diff)
            })
            .reduce(|| (0.0, 0.0), |(n1, d1), (n2, d2)| (n1 + n2, d1 + d2));

        return Some(numerator / denominator);
    }
}

#[pyclass]
pub struct PiecewisePolynomial {
    pub x_ranges: Vec<(f64, f64)>,
    pub y_functions: Vec<Polynomial>,
}

impl Evaluatable for PiecewisePolynomial {
    fn eval(&self, x: f64) -> Option<f64> {
        if x < self.x_ranges[0].0 {
            return self.y_functions[0].eval(x);
        }
        for i in 0..self.x_ranges.len() {
            let x_min = self.x_ranges[i].0;
            let x_max = self.x_ranges[i].1;
            if x >= x_min && x <= x_max {
                return self.y_functions[i].eval(x);
            }
        }
        return self.y_functions[self.y_functions.len() - 1].eval(x);
    }
}

impl Differentiable for PiecewisePolynomial {
    fn differentiate(&self) -> Self {
        let mut diff_poly = Vec::new();
        for polynomial in self.y_functions.iter() {
            diff_poly.push(polynomial.differentiate())
        }
        return PiecewisePolynomial {
            x_ranges: self.x_ranges.clone(),
            y_functions: diff_poly,
        };
    }
}

impl Integrable for PiecewisePolynomial {
    fn integrate(&self, x0: f64, y0: f64) -> PiecewisePolynomial {
        let (i, p) = if x0 < self.x_ranges[0].0 {
            (0, self.y_functions[0].integrate(x0, y0))
        } else if x0 > self.x_ranges[self.x_ranges.len() - 1].1 {
            (
                self.x_ranges.len() - 1,
                self.y_functions[self.y_functions.len() - 1].integrate(x0, y0),
            )
        } else {
            let mut i: usize = 0;
            while i < self.x_ranges.len() {
                let x_min = self.x_ranges[i].0;
                let x_max = self.x_ranges[i].1;
                if x0 >= x_min && x0 <= x_max {
                    break;
                };
                i += 1;
            }
            (i, self.y_functions[i].integrate(x0, y0))
        };
        let mut integrated_functions = self.y_functions.clone();
        integrated_functions[i] = p.clone();

        // integrate from i-1 to 0
        let mut j = (i - 1) as i32;
        let mut x_lower = self.x_ranges[i].0;
        let mut y_lower = p.eval(x_lower).unwrap();
        while j >= 0 {
            let p_lower = self.y_functions[j as usize].integrate(x_lower, y_lower);
            x_lower = self.x_ranges[j as usize].0;
            y_lower = p_lower.eval(x_lower).unwrap();
            integrated_functions[j as usize] = p_lower;
            j -= 1;
        }
        // integrate from i+1 to n
        let mut k = i + 1;
        let mut x_upper = self.x_ranges[i].1;
        let mut y_upper = p.eval(x_upper).unwrap();
        while k < self.x_ranges.len() {
            let p_upper = self.y_functions[k].integrate(x_upper, y_upper);
            x_upper = self.x_ranges[k].1;
            y_upper = p_upper.eval(x_upper).unwrap();
            integrated_functions[k] = p_upper;
            k += 1;
        }
        return PiecewisePolynomial {
            x_ranges: self.x_ranges.clone(),
            y_functions: integrated_functions,
        };
    }
}
