mod root_finding;

use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

use crate::root_finding::herons_method::herons_method;
use crate::root_finding::secant_method::secant_method;
use crate::root_finding::bisection_method::bisection_method;
use crate::root_finding::newton_raphson_method::newton_raphson_method;


type Function = Box<dyn Fn(f64) -> f64>;


fn wrap_py_function(py_function: Py<PyAny>) -> Function {
    // wrap python function to rust function on heap
    Box::new(move |x| {
        let y: f64 = Python::attach(
            |py| py_function.call1(py, (x,)).unwrap().as_ref().extract(py).unwrap());
        y
    })
}


#[pyfunction(name = "herons_method")]
#[pyo3(signature = (a, x_0 = 1.0, n_max=100))]
pub fn herons_method_py(a: f64, x_0: f64, n_max: i64) -> f64 {
    herons_method(a, x_0, n_max)
}


#[pyfunction(name = "bisection_method")]
#[pyo3(signature = (f, a, b, n_max=100, eps_tol=0.000001))]
pub fn bisection_method_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    n_max: i64,
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
    n_max: i64,
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
    n_max: i64,
    eps_tol: f64
) -> PyResult<f64> {
    let f: Function = wrap_py_function(function);
    Ok(secant_method(f, x_0, x_1, n_max, eps_tol))
}


#[pymodule]
fn comp_math(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(herons_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(secant_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(bisection_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(newton_raphson_method_py, m)?)?;
    Ok(())
}
