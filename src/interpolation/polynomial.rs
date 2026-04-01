use pyo3::prelude::*;


#[pyclass]
#[derive(Debug)]
pub struct Polynomial {
    coefficients: Vec<f64>,
}


pub(crate) trait Evaluatable {
    fn eval(&self, x: f64) -> f64;
}


trait PythonCallable {
    fn __call__(&self, x: f64) -> f64;
}


impl Polynomial {
    pub(crate) fn new(coefficients: Vec<f64>) -> Polynomial {
        Polynomial { coefficients }
    }
}


impl Evaluatable for Polynomial {
    fn eval(&self, x: f64) -> f64 {
        let mut x_power = 1.0;
        let mut total = 0.0;
        for coefficient in &self.coefficients {
            total += coefficient * x_power;
            x_power *= x;
        }
        total
    }
}


impl PythonCallable for Polynomial where Polynomial: Evaluatable {
    fn __call__(&self, x: f64) -> f64 {
        self.eval(x)
    }
}
