use pyo3::pyclass;

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


impl Differentiable for Polynomial {
    fn differentiate(&self) -> Self {
        let mut diff_weights = Vec::new();
        for i in 1..self.weights.len() {
            diff_weights.push(self.weights[i] * i as f64);
        }

        Polynomial{
            weights: diff_weights,
            x_i: self.x_i,
        }
    }
}


impl Integrable for Polynomial {
    fn integrate(&self, x0: f64, y0: f64) -> Self {
        let mut int_weights = Vec::new();
        int_weights.push(0.0);
        for i in 0..self.weights.len() {
            int_weights.push(self.weights[i] / (i+1) as f64);
        }

        let mut definite_integral = Polynomial{
            weights: int_weights,
            x_i: self.x_i,
        };
        let w_0 = y0 - definite_integral.eval(x0).unwrap();
        definite_integral.weights[0] = w_0;
        definite_integral
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
        Some(total)
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
        let mut numerator = 0.0;
        let mut denominator = 0.0;
        for i in 0..self.xs.len() {
            numerator += (self.weights[i] * self.ys[i])/(x - self.xs[i]);
            denominator += self.weights[i]/(x - self.xs[i]);
        }
        Some(numerator/denominator)
    }
}