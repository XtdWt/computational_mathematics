mod herons_method;

use pyo3::prelude::*;
use herons_method::herons_method;


#[pyfunction(name="herons_method")]
#[pyo3(signature = (a, x_0 = 1.0, n_max=100))]
pub fn herons_method_py(a: f64, x_0: f64, n_max: i64) -> f64 {
    herons_method(a, x_0, n_max)
}


#[pymodule]
fn comp_math(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(herons_method_py, m)?)?;
    Ok(())
}
