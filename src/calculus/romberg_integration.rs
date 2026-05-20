use crate::Function;
use nalgebra::DMatrix;

pub fn romberg_integration(
    f: &Function,
    a: f64,
    b: f64,
    k: usize,
) -> f64 {
    let interval_len = b - a;
    let mut r_matrix = DMatrix::<f64>::zeros(k, k);
    r_matrix[(0, 0)] = interval_len * (f(a) - f(b)) / 2.0;

    for j in 1..k {
        let h_j = interval_len / 2.0_f64.powi(j as i32);
        let a1 = r_matrix[(j-1, 0)] / 2.0;
        let end = 2_i32.pow((j as u32 ) - 1) + 1;
        let a2: f64 = (1..end as usize)
            .map(|i| f(a + (2.0 * (i as f64) - 1.0) * h_j))
            .sum();
        r_matrix[(j, 0)] = a1 + h_j * a2;

        for k in 1..j+1 {
            let numer = 4.0_f64.powi(k as i32) * r_matrix[(j, k-1)] - r_matrix[(j-1, k-1)];
            let denom = 4.0_f64.powi(k as i32) - 1.0;
            r_matrix[(j, k)] = numer / denom;
        }
    }
    return r_matrix[(k-1, k-1)];
}


#[cfg(test)]
mod tests {
    use super::*;
    use crate::assert_almost_eq;

    #[test]
    fn test_simple_polynomial_k_two() {
        let f: Function = Box::new(|x| x.powi(2) + 7.0);
        let r = romberg_integration(&f, -1.0, 4.0, 2);

        assert_almost_eq!(r, 55.0/3.0);
    }

    #[test]
    fn test_simple_polynomial_k_large() {
        let k = 25;
        let f: Function = Box::new(|x| x.powi(2) + 7.0);
        let r = romberg_integration(&f, -1.0, 4.0, k);

        assert_almost_eq!(r, 170.0 / 3.0, 1e-4);
    }
}
