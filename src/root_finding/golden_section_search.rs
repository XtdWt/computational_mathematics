use std::f64::consts::GOLDEN_RATIO;


pub fn golden_section_search<F: Fn(f64) -> f64>(
    f: &F,
    a: f64,
    b: f64,
    n_max: usize,
    eps_tol: f64,
) -> f64 {
    let mut a_i = a;
    let mut b_i = b;
    for _ in 0..n_max {
        if b_i - a_i < eps_tol {
            return (b_i + a_i) / 2.0
        }

        let h = b_i - a_i;
        let x_1 = b_i - h * GOLDEN_RATIO.recip();
        let x_2 = a_i + h * GOLDEN_RATIO.recip();
        if f(x_1) < f(x_2) {
            b_i = x_2;
        } else {
            a_i = x_1;
        }
    }
    return (b_i + a_i) / 2.0
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;
    use std::f64::consts::PI;

    #[test]
    fn test_simple_golden_section_search() {
        let f = |x: f64| x.powi(2) - 1.0;
        let x_star = golden_section_search(&f, -2.0, 1.0, 20, 1e-6);
        assert_almost_eq!(x_star, 0.0, 1e-3)
    }

    #[test]
    fn test_simple_golden_section_search_early_stop() {
        let f = |x: f64| x.powi(2) - 1.0;
        let x_star = golden_section_search(&f, -3.0, 3.0, 1_000, 1e-100);
        assert_almost_eq!(x_star, 0.0)
    }

    #[test]
    fn test_sine_golden_section_search() {
        let f = |x: f64| x.cos();
        let x_star = golden_section_search(&f, 0.0, PI, 20, 1e-8);
        assert_almost_eq!(x_star, PI, 1e-3)
    }
}
