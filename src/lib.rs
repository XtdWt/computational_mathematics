use num_complex::Complex;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

mod root_finding;
use crate::root_finding::herons_method::herons_method;
use crate::root_finding::secant_method::secant_method;
use crate::root_finding::bisection_method::bisection_method;
use crate::root_finding::newton_raphson_method::newton_raphson_method;

mod interpolation;
use crate::interpolation::polynomial::{
    Polynomial,
    PiecewisePolynomial,
    LagrangePolynomial,
    NewtonsDividedDifferencePolynomial,
    Differentiable,
    Integrable,
    Evaluatable,
};
use crate::interpolation::barycentric_lagrange_interpolation::{
    barycentric_lagrange_interpolation,
};
use crate::interpolation::newtons_divided_difference_interpolation::{
    newtons_divided_difference_interpolation
};
use crate::interpolation::chebyshev_nodes::chebyshev_nodes;
use crate::interpolation::cubic_spline_interpolation::{
    cubic_spline_interpolation,
};
use crate::interpolation::fast_fourier_transform::{
    fast_fourier_transform,
    inverse_fast_fourier_transform,
    fast_fourier_transform_frequencies,
};

type Function = Box<dyn Fn(f64) -> f64>;


fn wrap_py_function(py_function: Py<PyAny>) -> Function {
    // wrap python function to rust function on heap
    Box::new(move |x| {
        let y: f64 = Python::attach(
            |py| py_function.call1(py, (x,)).unwrap().extract(py).unwrap());
        y
    })
}


#[pyfunction(name = "herons_method")]
#[pyo3(signature = (a, x_0 = 1.0, n_max=100))]
pub fn herons_method_py(a: f64, x_0: f64, n_max: usize) -> f64 {
    herons_method(a, x_0, n_max)
}


#[pyfunction(name = "bisection_method")]
#[pyo3(signature = (f, a, b, n_max=100, eps_tol=0.000001))]
pub fn bisection_method_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f: Function = wrap_py_function(f);
    if f(a) * f(b) > 0.0 {
        return Err(PyValueError::new_err(
            "`f(a)` and `f(b)` have the same sign.",
        ));
    }
    Ok(bisection_method(f, a, b, n_max, eps_tol))
}


#[pyfunction(name = "newton_raphson_method")]
#[pyo3(signature = (function, derivative, x_0, n_max=100, eps_tol=0.000001))]
pub fn newton_raphson_method_py(
    function: Py<PyAny>,
    derivative: Py<PyAny>,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f: Function = wrap_py_function(function);
    let df: Function = wrap_py_function(derivative);
    Ok(newton_raphson_method(f, df, x_0, n_max, eps_tol))
}


#[pyfunction(name = "secant_method")]
#[pyo3(signature = (function, x_0, x_1, n_max=100, eps_tol=0.000001))]
pub fn secant_method_py(
    function: Py<PyAny>,
    x_0: f64,
    x_1: f64,
    n_max: usize,
    eps_tol: f64
) -> PyResult<f64> {
    let f: Function = wrap_py_function(function);
    Ok(secant_method(f, x_0, x_1, n_max, eps_tol))
}


#[pymethods]
impl LagrangePolynomial {
    fn __call__(&self, x:f64) -> f64 {
        self.eval(x).unwrap_or_else(|| f64::NAN)
    }
}


#[pyfunction(name = "barycentric_lagrange_interpolation")]
pub fn barycentric_lagrange_interpolation_py(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> LagrangePolynomial {
    barycentric_lagrange_interpolation(xs, ys)
}


#[pymethods]
impl NewtonsDividedDifferencePolynomial {
    fn __call__(&self, x:f64) -> f64 {
        self.eval(x).unwrap_or_else(|| f64::NAN)
    }
}


#[pyfunction(name = "newtons_divided_difference_interpolation")]
pub fn newtons_divided_difference_interpolation_py(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> NewtonsDividedDifferencePolynomial {
    newtons_divided_difference_interpolation(xs, ys)
}


#[pyfunction(name = "chebyshev_nodes")]
pub fn chebyshev_nodes_py(
    a: f64,
    b: f64,
    n: usize,
) -> Vec<f64> {
    chebyshev_nodes(a, b, n)
}

#[pymethods]
impl Polynomial {
    fn __call__(&self, x: f64) -> f64 {
        self.eval(x).unwrap_or_else(|| f64::NAN)
    }

    #[pyo3(name = "differentiate")]
    fn differentiate_py(&self) -> Self {
        self.differentiate()
    }

    #[pyo3(name = "integrate")]
    fn integrate_py(&self, x0: f64, y0: f64) -> Self {
        self.integrate(x0, y0)
    }
}


#[pymethods]
impl PiecewisePolynomial {
    fn __call__(&self, x: f64) -> f64 {
        self.eval(x).unwrap_or_else(|| f64::NAN)
    }

    #[pyo3(name = "differentiate")]
    fn differentiate_py(&self) -> Self {
        self.differentiate()
    }

    #[pyo3(name = "integrate")]
    fn integrate_py(&self, x0: f64, y0: f64) -> Self {
        self.integrate(x0, y0)
    }
}


#[pyfunction(name = "cubic_spline_interpolation")]
pub fn cubic_spline_interpolation_py(xs: Vec<f64>, ys: Vec<f64>) -> PiecewisePolynomial {
    cubic_spline_interpolation(xs, ys)
}


#[pyfunction(name = "fast_fourier_transform")]
pub fn fast_fourier_transform_py(
    xs: Vec<Complex<f64>>,
) -> Vec<Complex<f64>> {
    fast_fourier_transform(xs)
}


#[pyfunction(name = "inverse_fast_fourier_transform")]
pub fn inverse_fast_fourier_transform_py(
    xs: Vec<Complex<f64>>,
) -> Vec<Complex<f64>> {
    inverse_fast_fourier_transform(xs)
}


#[pyfunction(name = "fast_fourier_transform_frequencies")]
pub fn fast_fourier_transform_frequencies_py(
    n: usize,
    d: f64,
) -> Vec<f64> {
    fast_fourier_transform_frequencies(n, d)
}


#[pymodule]
fn computational_mathematics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(herons_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(secant_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(bisection_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(newton_raphson_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(barycentric_lagrange_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(newtons_divided_difference_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(chebyshev_nodes_py, m)?)?;
    m.add_function(wrap_pyfunction!(cubic_spline_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(fast_fourier_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(inverse_fast_fourier_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fast_fourier_transform_frequencies_py, m)?)?;

    m.add_class::<PiecewisePolynomial>()?;
    m.add_class::<Polynomial>()?;
    m.add_class::<LagrangePolynomial>()?;
    m.add_class::<NewtonsDividedDifferencePolynomial>()?;

    Ok(())
}
