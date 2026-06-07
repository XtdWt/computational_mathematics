use num_complex::Complex;

use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

mod root_finding;
use crate::root_finding::bisection_method::bisection_method;
use crate::root_finding::herons_method::herons_method;
use crate::root_finding::newton_raphson_method::{newton_raphson_method, quasi_newton_raphson_method};
use crate::root_finding::secant_method::secant_method;
use crate::root_finding::golden_section_search::golden_section_search;

mod interpolation;
use crate::interpolation::barycentric_lagrange_interpolation::barycentric_lagrange_interpolation;
use crate::interpolation::chebyshev_nodes::chebyshev_nodes;
use crate::interpolation::cubic_spline_interpolation::cubic_spline_interpolation;
use crate::interpolation::fast_fourier_transform::{
    fast_fourier_transform, fast_fourier_transform_frequencies, inverse_fast_fourier_transform,
};
use crate::interpolation::newtons_divided_difference_interpolation::newtons_divided_difference_interpolation;
use crate::interpolation::polynomial::{
    Differentiable, Evaluatable, Integrable, LagrangePolynomial,
    NewtonsDividedDifferencePolynomial, PiecewisePolynomial, Polynomial,
};

mod calculus;
use crate::calculus::first_derivative::first_derivative;
use crate::calculus::composite_integration::{composite_simpsons_rule, composite_trapezoid_rule};
use crate::calculus::second_derivative::second_derivative;
use crate::calculus::adaptive_quadrature::adaptive_quadrature;
use crate::calculus::romberg_integration::romberg_integration;
use crate::calculus::simple_ivp_solver::simple_ivp_solver;
use crate::calculus::rk4_solver::rk4_solver;
use crate::calculus::util::{DerivativeType, SimpleIVPSolverType};

pub fn wrap_py_multivariatefunction(py_function: Py<PyAny>) -> impl Fn(f64, Vec<f64>) -> Vec<f64> {
    move |t: f64, x: Vec<f64>| {
        Python::attach(|py| {
            py_function
                .call1(py, (t, x,))
                .and_then(|result| result.extract::<Vec<f64>>(py))
                .unwrap_or_else(|e| {
                    e.restore(py);
                    return Vec::new();
                })
        })
    }
}

pub fn wrap_py_function(py_function: Py<PyAny>) -> impl Fn(f64) -> f64 {
    move |x: f64| {
        Python::attach(|py| {
            py_function
                .call1(py, (x,))
                .and_then(|result| result.extract::<f64>(py))
                .unwrap_or_else(|e| {
                    e.restore(py);
                    return f64::NAN;
                })
        })
    }
}

#[macro_export]
macro_rules! assert_almost_eq {
    ($a:expr, $b:expr, $eps:expr) => {
        if ($a - $b).abs() > $eps {
            panic!(
                "assertion `|left - right| <= eps` failed for\n left={} \n right={} \n difference={} > epsilon={}", $a, $b, ($a - $b).abs(), $eps
            )
        }
    };

    //epsilon default value of 1e-7
    ($a:expr, $b:expr) => {
        if ($a - $b).abs() > 1e-7 {
            panic!(
                "assertion `|left - right| <= eps` failed for\n left={} \n right={} \n difference={} > epsilon={}", $a, $b, ($a - $b).abs(), 1e-7
            )
        }
    };
}


#[pyfunction(name = "herons_method")]
#[pyo3(signature = (a, x_0 = 1.0, n_max=100))]
pub fn herons_method_py(a: f64, x_0: f64, n_max: usize) -> f64 {
    return herons_method(a, x_0, n_max);
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
    let f = wrap_py_function(f);
    let f_a = f(a);
    let f_b = f(b);
    if f_a * f_b > 0.0 {
        return Err(PyValueError::new_err(
            "`f(a)` and `f(b)` have the same sign.",
        ));
    }
    return Ok(bisection_method(&f, a, b, n_max, eps_tol));
}

#[pyfunction(name = "newton_raphson_method")]
#[pyo3(signature = (f, df, x_0, n_max=100, eps_tol=0.000001))]
pub fn newton_raphson_method_py(
    f: Py<PyAny>,
    df: Py<PyAny>,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    let df = wrap_py_function(df);
    return Ok(newton_raphson_method(&f, &df, x_0, n_max, eps_tol));
}

#[pyfunction(name = "quasi_newton_raphson_method")]
#[pyo3(signature = (f, x_0, n_max=100, eps_tol=0.000001))]
pub fn quasi_newton_raphson_method_py(
    f: Py<PyAny>,
    x_0: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(quasi_newton_raphson_method(&f, x_0, n_max, eps_tol));
}

#[pyfunction(name = "secant_method")]
#[pyo3(signature = (f, x_0, x_1, n_max=100, eps_tol=0.000001))]
pub fn secant_method_py(
    f: Py<PyAny>,
    x_0: f64,
    x_1: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(secant_method(&f, x_0, x_1, n_max, eps_tol));
}

#[pyfunction(name = "golden_section_search")]
#[pyo3(signature = (f, a, b, n_max=100, eps_tol=0.000001))]
pub fn golden_section_search_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    n_max: usize,
    eps_tol: f64,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(golden_section_search(&f, a, b, n_max, eps_tol));
}

#[pymethods]
impl LagrangePolynomial {
    fn __call__(&self, x: f64) -> f64 {
        return self.eval(x).unwrap_or_else(|| f64::NAN);
    }
}

#[pyfunction(name = "barycentric_lagrange_interpolation")]
pub fn barycentric_lagrange_interpolation_py(xs: Vec<f64>, ys: Vec<f64>) -> LagrangePolynomial {
    return barycentric_lagrange_interpolation(xs, ys);
}

#[pymethods]
impl NewtonsDividedDifferencePolynomial {
    fn __call__(&self, x: f64) -> f64 {
        return self.eval(x).unwrap_or_else(|| f64::NAN);
    }
}

#[pyfunction(name = "newtons_divided_difference_interpolation")]
pub fn newtons_divided_difference_interpolation_py(
    xs: Vec<f64>,
    ys: Vec<f64>,
) -> NewtonsDividedDifferencePolynomial {
    return newtons_divided_difference_interpolation(xs, ys);
}

#[pyfunction(name = "chebyshev_nodes")]
pub fn chebyshev_nodes_py(a: f64, b: f64, n: usize) -> Vec<f64> {
    return chebyshev_nodes(a, b, n);
}

#[pymethods]
impl Polynomial {
    fn __call__(&self, py: Python<'_>, x: f64) -> f64 {
        return py.detach(|| { self.eval(x).unwrap_or_else(|| f64::NAN) });
    }

    #[pyo3(name = "differentiate")]
    fn differentiate_py(&self) -> Self {
        return self.differentiate();
    }

    #[pyo3(name = "integrate")]
    fn integrate_py(&self, x0: f64, y0: f64) -> Self {
        return self.integrate(x0, y0);
    }
}

#[pymethods]
impl PiecewisePolynomial {
    fn __call__(&self, x: f64) -> f64 {
        return self.eval(x).unwrap_or_else(|| f64::NAN);
    }

    #[pyo3(name = "differentiate")]
    fn differentiate_py(&self) -> Self {
        return self.differentiate();
    }

    #[pyo3(name = "integrate")]
    fn integrate_py(&self, x0: f64, y0: f64) -> Self {
        return self.integrate(x0, y0);
    }
}

#[pyfunction(name = "cubic_spline_interpolation")]
pub fn cubic_spline_interpolation_py(xs: Vec<f64>, ys: Vec<f64>) -> PiecewisePolynomial {
    return cubic_spline_interpolation(xs, ys);
}

#[pyfunction(name = "fast_fourier_transform")]
pub fn fast_fourier_transform_py(py: Python<'_>, xs: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    return py.detach(move || fast_fourier_transform(xs));
}

#[pyfunction(name = "inverse_fast_fourier_transform")]
pub fn inverse_fast_fourier_transform_py(xs: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    return inverse_fast_fourier_transform(xs);
}

#[pyfunction(name = "fast_fourier_transform_frequencies")]
pub fn fast_fourier_transform_frequencies_py(n: usize, d: f64) -> Vec<f64> {
    return fast_fourier_transform_frequencies(n, d);
}

#[pyfunction(name = "first_derivative")]
#[pyo3(signature = (f, x, h, method="central"))]
pub fn first_derivative_py(f: Py<PyAny>, x: f64, h: f64, method: &str) -> PyResult<f64> {
    let f = wrap_py_function(f);
    let method_as_enum = if method == "backward" {
        DerivativeType::Backward
    } else if method == "forward" {
        DerivativeType::Forward
    } else if method == "central" {
        DerivativeType::Central
    } else {
        return Err(PyValueError::new_err(
            "method: ".to_owned() + method + " must be one of ['backward', 'forward', 'central']",
        ));
    };
    return Ok(first_derivative(&f, x, h, method_as_enum));
}

#[pyfunction(name = "second_derivative")]
#[pyo3(signature = (f, x, h, method="central"))]
pub fn second_derivative_py(f: Py<PyAny>, x: f64, h: f64, method: &str) -> PyResult<f64> {
    let f = wrap_py_function(f);
    let method_as_enum = if method == "backward" {
        DerivativeType::Backward
    } else if method == "forward" {
        DerivativeType::Forward
    } else if method == "central" {
        DerivativeType::Central
    } else {
        return Err(PyValueError::new_err(
            "method: ".to_owned() + method + " must be one of ['backward', 'forward', 'central']",
        ));
    };
    return Ok(second_derivative(&f, x, h, method_as_enum));
}

#[pyfunction(name = "composite_trapezoid_rule")]
#[pyo3(signature = (f, a, b, n_buckets=100))]
pub fn composite_trapezoid_rule_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(composite_trapezoid_rule(&f, a, b, n_buckets));
}

#[pyfunction(name = "composite_simpsons_rule")]
#[pyo3(signature = (f, a, b, n_buckets=100))]
pub fn composite_simpsons_rule_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    n_buckets: usize,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(composite_simpsons_rule(&f, a, b, n_buckets));
}

#[pyfunction(name = "adaptive_quadrature")]
#[pyo3(signature = (f, a, b, eps_tol=1e-7))]
pub fn adaptive_quadrature_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    eps_tol: f64,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(adaptive_quadrature(&f, a, b, eps_tol));
}


#[pyfunction(name = "romberg_integration")]
#[pyo3(signature = (f, a, b, k))]
pub fn romberg_integration_py(
    f: Py<PyAny>,
    a: f64,
    b: f64,
    k: usize,
) -> PyResult<f64> {
    let f = wrap_py_function(f);
    return Ok(romberg_integration(&f, a, b, k));
}


#[pyfunction(name = "simple_ivp_solver")]
#[pyo3(signature = (df, y0, t0, tn, h=1e-6, method="midpoint"))]
pub fn simple_ivp_solver_py(
    df: Py<PyAny>,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
    method: &str,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let df = wrap_py_multivariatefunction(df);
    let method_as_enum = if method == "euler" {
        SimpleIVPSolverType::Eulers
    } else if method == "trapeziod" {
        SimpleIVPSolverType::Trapeziod
    } else if method == "midpoint" {
        SimpleIVPSolverType::Midpoint
    } else {
        return Err(PyValueError::new_err(
            "method: ".to_owned() + method + " must be one of ['euler', 'trapeziod', 'midpoint']",
        ));
    };
    return Ok(simple_ivp_solver(&df, y0, t0, tn, h, method_as_enum));
}

#[pyfunction(name = "rk4_solver")]
#[pyo3(signature = (df, y0, t0, tn, h=1e-6))]
pub fn rk4_solver_py(
    df: Py<PyAny>,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
) -> PyResult<(Vec<f64>, Vec<Vec<f64>>)> {
    let df = wrap_py_multivariatefunction(df);
    return Ok(rk4_solver(&df, y0, t0, tn, h));
}

#[pymodule]
fn computational_mathematics(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(herons_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(secant_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(bisection_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(newton_raphson_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(quasi_newton_raphson_method_py, m)?)?;
    m.add_function(wrap_pyfunction!(golden_section_search_py, m)?)?;

    m.add_function(wrap_pyfunction!(barycentric_lagrange_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(newtons_divided_difference_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(chebyshev_nodes_py, m)?)?;
    m.add_function(wrap_pyfunction!(cubic_spline_interpolation_py, m)?)?;
    m.add_function(wrap_pyfunction!(fast_fourier_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(inverse_fast_fourier_transform_py, m)?)?;
    m.add_function(wrap_pyfunction!(fast_fourier_transform_frequencies_py, m)?)?;

    m.add_function(wrap_pyfunction!(first_derivative_py, m)?)?;
    m.add_function(wrap_pyfunction!(second_derivative_py, m)?)?;
    m.add_function(wrap_pyfunction!(composite_trapezoid_rule_py, m)?)?;
    m.add_function(wrap_pyfunction!(composite_simpsons_rule_py, m)?)?;
    m.add_function(wrap_pyfunction!(adaptive_quadrature_py, m)?)?;
    m.add_function(wrap_pyfunction!(romberg_integration_py, m)?)?;
    m.add_function(wrap_pyfunction!(simple_ivp_solver_py, m)?)?;
    m.add_function(wrap_pyfunction!(rk4_solver_py, m)?)?;

    m.add_class::<PiecewisePolynomial>()?;
    m.add_class::<Polynomial>()?;
    m.add_class::<LagrangePolynomial>()?;
    m.add_class::<NewtonsDividedDifferencePolynomial>()?;

    return Ok(());
}
