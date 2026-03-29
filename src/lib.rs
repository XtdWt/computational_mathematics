mod root_finding;

use pyo3::prelude::*;
use crate::root_finding::herons_method::herons_method;
use crate::root_finding::bisection_method::bisection_method;


type F = Box<dyn Fn(f64) -> f64>;


#[pyfunction(name="herons_method")]
#[pyo3(signature = (a, x_0 = 1.0, n_max=100))]
pub fn herons_method_py(a: f64, x_0: f64, n_max: i64) -> f64 {
    herons_method(a, x_0, n_max)
}


#[pyfunction(name="bisection_method")]
#[pyo3(signature = (function, a, b, n_max=100))]
pub fn bisection_method_py(function: Py<PyAny>, a: f64, b: f64, n_max: i64) -> Option<f64> {
    // wrap python function to rust function on heap
    let rust_function: F = Box::new(move |x| {
        Python::attach(|py| {
            function.call1(py, (x,)).unwrap();
        });
    x});
    bisection_method(rust_function, a, b, n_max)
}


#[pymodule]
fn comp_math(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(herons_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(bisection_method_py, m)?)?;
    Ok(())
}
