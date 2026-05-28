use crate::MultivariateFunction;
use crate::calculus::util::SimpleIVPSolverType;

fn t_range(t_start: f64, t_end: f64, step: f64) -> Vec<f64> {
    let size = (( t_end - t_start ) / step) as usize;
    return (0..size + 1)
        .map(|i| t_start + i as f64 * step)
        .collect();
}

fn eulers_ivp(
    df: &MultivariateFunction,
    y0: f64,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<f64> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        // only take the second `y` value from df(t, y)
        ys.push(ys[i-1] + h * df( vec![ts[i-1], ys[i-1]] )[0]);
    }
    return (ts, ys);
}

fn trapeziod_ivp(
    df: &MultivariateFunction,
    y0: f64,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<f64> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        // only take the second `y` value from df(t, y)
        let df_ti = df( vec![ts[i-1], ys[i-1]] )[0];
        ys.push(ys[i-1] + (h/2.0)*(df_ti + df( vec![ts[i], ys[i-1] + h*df_ti] )[0]));
    }
    return (ts, ys);
}

fn midpoint_ivp(
    df: &MultivariateFunction,
    y0: f64,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<f64>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<f64> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        // only take the second `y` value from df(t, y)
        let df_ti = df( vec![ts[i-1], ys[i-1]] )[0];
        ys.push(ys[i-1] + h * df( vec![ts[i-1] + h/2.0, ys[i-1] + (h/2.0) * df_ti] )[0]);
    }
    return (ts, ys);
}

pub fn simple_ivp_solver(
    df: &MultivariateFunction,
    y0: f64,
    t0: f64,
    tn: f64,
    h: f64,
    method: SimpleIVPSolverType,
) -> (Vec<f64>, Vec<f64>) {
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
        let df: MultivariateFunction = Box::new(|x| vec![0.5 + 0.5*x[1].powi(2)]);

        let (ts, ys) = simple_ivp_solver(&df, 3_f64.sqrt(), 0.0, 1.0, 1e-6, SimpleIVPSolverType::Eulers);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1], (0.5 + PI/3.0).tan(), 1e-2);
        let (ts, ys) = simple_ivp_solver(&df, 3_f64.sqrt(), 0.0, 1.0, 1e-6, SimpleIVPSolverType::Trapeziod);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1], (0.5 + PI/3.0).tan(), 1e-2);
        let (ts, ys) = simple_ivp_solver(&df, 3_f64.sqrt(), 0.0, 1.0, 1e-6, SimpleIVPSolverType::Midpoint);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1], (0.5 + PI/3.0).tan(), 1e-2);
    }
}
