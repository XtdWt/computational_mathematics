pub fn collocation_method<
    F: Fn(f64) -> f64,
    DF: Fn(f64, Vec<f64>) -> Vec<f64>,
>(
    df: &DF,
    boundary: (f64, f64),

) -> {}
