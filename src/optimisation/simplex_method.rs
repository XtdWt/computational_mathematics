use nalgebra::DMatrix;

pub fn simplex_method(
    c_values: Vec<f64>,
    a_values: Vec<Vec<f64>>,
    b_values: Vec<f64>,
    x_bounds: Vec<f64>,
) -> (f64, Vec<f64>) {
    let c = DMatrix::from_vec(c_values.len(), 1, c_values);
    let a = DMatrix::from_row_slice(
        b_values.len(), a_values.len(), &(a_values
            .into_iter()
            .flatten()
            .collect::<Vec<f64>>()
        )
    );

    let b = DMatrix::from_vec(b_values.len(), 1, b_values);

    let x = a.try_inverse().unwrap() * b;

    let cost = c.dot(&x);


    return (cost, x.as_slice().to_vec());
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_simple() {
        let c_values = vec![1.0, 2.0, 3.0];
        let a_values = vec![
            vec![1.0, 0.0, 0.0],
            vec![0.0, 1.0, 0.0],
            vec![0.0, 0.0, 1.0],
        ];
        let b_values = vec![1.0, 2.0, 3.0];

        let (cost, x) = simplex_method(
            c_values, a_values, b_values, vec![]
        );

        assert_eq!(cost, 14.0);
        assert_eq!(x, vec![1.0, 2.0, 3.0]);
    }
}
