use num_complex::Complex;
use std::f64::consts::PI;

fn fast_fourier_transform_base(xs: Vec<Complex<f64>>, is_inverse: bool) -> Vec<Complex<f64>> {
    let n = xs.len();
    if n <= 1 {
        return xs;
    };
    let x_hat_even =
        fast_fourier_transform_base(xs.iter().step_by(2).cloned().collect(), is_inverse);
    let x_hat_odd =
        fast_fourier_transform_base(xs.iter().skip(1).step_by(2).cloned().collect(), is_inverse);

    let mut x_hat: Vec<Complex<f64>> = vec![Complex::from(0.0); n];
    for i in 0..n / 2 {
        let w = if is_inverse {
            Complex::from_polar(1.0, (2.0 * PI * (i as f64)) / (n as f64))
        } else {
            Complex::from_polar(1.0, (-2.0 * PI * (i as f64)) / (n as f64))
        };
        x_hat[i] = x_hat_even[i] + w * x_hat_odd[i];
        x_hat[i + n / 2] = x_hat_even[i] - w * x_hat_odd[i];
    }
    return x_hat;
}

fn pad_zeroes_complex(mut xs: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = xs.len();
    if !(n > 0 && n.count_ones() == 1) {
        let next_pow_two = n.next_power_of_two();
        for _ in n..next_pow_two {
            xs.push(Complex::new(0.0, 0.0));
        }
    }
    xs
}

pub fn fast_fourier_transform(xs: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let xs = pad_zeroes_complex(xs);
    fast_fourier_transform_base(xs, false)
}

pub fn inverse_fast_fourier_transform(xs: Vec<Complex<f64>>) -> Vec<Complex<f64>> {
    let n = xs.len() as f64;
    let xs = pad_zeroes_complex(xs);
    fast_fourier_transform_base(xs, true)
        .iter()
        .map(|&x| x / n)
        .collect()
}

pub fn fast_fourier_transform_frequencies(n: usize, d: f64) -> Vec<f64> {
    let divisor = (n as f64) * d;
    let mut ts = Vec::new();
    let mid_value = if n % 2 == 0 { n / 2 } else { n / 2 + 1 };

    for i in 0..mid_value {
        ts.push((i as f64) / divisor);
    }
    for i in (1..(n / 2) + 1).rev() {
        ts.push(-(i as f64) / divisor);
    }
    ts
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fast_fourier_transform_four_points() {
        let test_xs = vec![
            Complex::from(2.0),
            Complex::from(3.0),
            Complex::from(2.0),
            Complex::from(3.0),
        ];

        let result = fast_fourier_transform(test_xs.clone());
        assert_eq!(
            result,
            vec![
                Complex::new(10.0, 0.0),
                Complex::new(0.0, 0.0),
                Complex::new(-2.0, 0.0),
                Complex::new(0.0, 0.0),
            ]
        );

        let inverse_result = inverse_fast_fourier_transform(result);
        assert_eq!(inverse_result, test_xs);
    }

    #[test]
    fn test_fast_fourier_transform_eight_points() {
        let test_xs = vec![
            Complex::from(2.0),
            Complex::from(3.0),
            Complex::from(2.0),
            Complex::from(3.0),
            Complex::from(4.0),
            Complex::from(5.0),
            Complex::from(0.0),
            Complex::from(0.0),
        ];

        let result = fast_fourier_transform(test_xs.clone());

        let inverse_result = inverse_fast_fourier_transform(result);
        let almost_inverse_result: Vec<Complex<f64>> = inverse_result
            .iter()
            // set "small" values to zero for ease of testing
            .map(|&x| {
                let re_part = if x.re.abs() < 1e-8 { 0.0 } else { x.re };
                let im_part = if x.im.abs() < 1e-8 { 0.0 } else { x.im };
                Complex::new(re_part, im_part)
            })
            .collect();
        assert_eq!(almost_inverse_result, test_xs);
    }

    #[test]
    fn test_padding_works() {
        let mut test_xs = vec![
            Complex::from(2.0),
            Complex::from(3.0),
            Complex::from(2.0),
            Complex::from(3.0),
            Complex::from(4.0),
            Complex::from(5.0),
        ];
        let result_six_points = fast_fourier_transform(test_xs.clone());

        test_xs.push(Complex::from(0.0));
        test_xs.push(Complex::from(0.0));
        let result_eight_points = fast_fourier_transform(test_xs.clone());

        assert_eq!(result_six_points, result_eight_points)
    }

    #[test]
    fn test_fast_fourier_transform_frequencies_odd() {
        let result = fast_fourier_transform_frequencies(5, 1.0);
        assert_eq!(result, vec![0.0, 0.2, 0.4, -0.4, -0.2]);

        let result = fast_fourier_transform_frequencies(3, 3.0_f64.recip());
        assert_eq!(result, vec![0.0, 1.0, -1.0]);

        let result = fast_fourier_transform_frequencies(1, 1.0);
        assert_eq!(result, vec![0.0]);

        let result = fast_fourier_transform_frequencies(7, 7.0_f64.recip());
        assert_eq!(result, vec![0.0, 1.0, 2.0, 3.0, -3.0, -2.0, -1.0]);
    }

    #[test]
    fn test_fast_fourier_transform_frequencies_even() {
        let result = fast_fourier_transform_frequencies(4, 1.0);
        assert_eq!(result, vec![0.0, 0.25, -0.5, -0.25]);

        let result = fast_fourier_transform_frequencies(2, 1.0);
        assert_eq!(result, vec![0.0, -0.5]);

        let result = fast_fourier_transform_frequencies(6, 3.0_f64.recip());
        assert_eq!(result, vec![0.0, 0.5, 1.0, -1.5, -1.0, -0.5]);

        let result = fast_fourier_transform_frequencies(8, 1.0);
        assert_eq!(
            result,
            vec![0.0, 0.125, 0.25, 0.375, -0.5, -0.375, -0.25, -0.125]
        );
    }
}
