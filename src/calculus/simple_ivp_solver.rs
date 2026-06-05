use crate::calculus::util::{SimpleIVPSolverType, t_range};

fn eulers_ivp<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: &DF,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<Vec<f64>> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        let delta_f = df(ts[i-1], ys[i-1].clone())
            .into_iter()
            .map(|x| h*x);
        ys.push(
            ys[i-1]
                .clone()
                .into_iter()
                .zip(delta_f)
                .map(|(a, b)| a + b)
                .collect()
        );
    }
    return (ts, ys);
}

fn trapeziod_ivp<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: &DF,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<Vec<f64>> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        let delta_f = df(ts[i-1], ys[i-1].clone());
        let delta_f_ti = df(
            ts[i],
            ys[i-1]
                .clone()
                .into_iter()
                .zip(delta_f.clone())
                .map(|(a, b)| a + h*b)
                .collect()
        );
        ys.push(
            ys[i-1]
                .clone()
                .into_iter()
                .zip(delta_f)
                .zip(delta_f_ti)
                .map(|((a, b), c)| a + (h/2.0) * (b + c))
                .collect()
        );
    }
    return (ts, ys);
}

fn midpoint_ivp<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: &DF,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<Vec<f64>> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        let delta_f = df(ts[i-1], ys[i-1].clone());
        ys.push(
            ys[i-1]
                .clone()
                .into_iter()
                .zip(
                    df(ts[i-1] + h/2.0, ys[i-1].clone().into_iter().zip(delta_f).map(|(a, b)| a + (h/2.0)*b).collect())
                )
                .map(|(a, b)| a + h*b)
                .collect()
        );
    }
    return (ts, ys);
}

pub fn simple_ivp_solver<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: &DF,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
    method: SimpleIVPSolverType,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    return match method {
        SimpleIVPSolverType::Eulers => eulers_ivp(&df, y0, t0, tn, h),
        SimpleIVPSolverType::Trapeziod => trapeziod_ivp(&df, y0, t0, tn, h),
        SimpleIVPSolverType::Midpoint => midpoint_ivp(&df, y0, t0, tn, h),
    };
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use crate::assert_almost_eq;

    #[test]
    fn test_simple_ivp_solver() {
        let df = |_t: f64, x: Vec<f64>| {vec![0.5 + 0.5*x[0].powi(2)]};

        let (ts, ys) = simple_ivp_solver(&df, vec![3_f64.sqrt()], 0.0, 1.0, 1e-6, SimpleIVPSolverType::Eulers);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1][0], (0.5 + PI/3.0).tan(), 1e-2);
        let (ts, ys) = simple_ivp_solver(&df, vec![3_f64.sqrt()], 0.0, 1.0, 1e-6, SimpleIVPSolverType::Trapeziod);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1][0], (0.5 + PI/3.0).tan(), 1e-2);
        let (ts, ys) = simple_ivp_solver(&df, vec![3_f64.sqrt()], 0.0, 1.0, 1e-6, SimpleIVPSolverType::Midpoint);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1][0], (0.5 + PI/3.0).tan(), 1e-2);
    }
}
