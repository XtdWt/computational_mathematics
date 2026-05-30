pub enum DerivativeType {
    Forward,
    Backward,
    Central,
}

pub enum SimpleIVPSolverType {
    Eulers,
    Trapeziod,
    Midpoint,
}

pub fn t_range(t_start: f64, t_end: f64, step: f64) -> Vec<f64> {
    let size = (( t_end - t_start ) / step) as usize;
    return (0..size + 1)
        .map(|i| t_start + i as f64 * step)
        .collect();
}
