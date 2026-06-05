use crate::calculus::util::t_range;


pub fn rk4_solver<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: DF,
    y0: Vec<f64>,
    t0: f64,
    tn: f64,
    h: f64,
) -> (Vec<f64>, Vec<Vec<f64>>) {
    let ts = t_range(t0, tn, h);
    let mut ys: Vec<Vec<f64>> = Vec::with_capacity(ts.len());
    ys.push(y0);
    for i in 1..ts.len() {
        let s1 = df(ts[i-1], ys[i-1].clone());
        let s2 = df(ts[i-1] + h/2.0, ys[i-1].clone().into_iter().zip(s1.clone()).map(|(a, b)| a + (h/2.0)*b).collect());
        let s3 = df(ts[i-1] + h/2.0, ys[i-1].clone().into_iter().zip(s2.clone()).map(|(a, b)| a + (h/2.0)*b).collect());
        let s4 = df(ts[i-1] + h, ys[i-1].clone().into_iter().zip(s3.clone()).map(|(a, b)| a + h*b).collect());

        ys.push(
            ys[i-1]
                .clone()
                .into_iter()
                .zip(s1)
                .zip(s2)
                .zip(s3)
                .zip(s4)
                .map(|((((y, s_1), s_2), s_3), s_4)| y + (h/6.0)*(s_1 + 2.0*s_2 + 2.0*s_3 + s_4))
                .collect()
        )
    }
    return (ts, ys);
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::PI;
    use crate::assert_almost_eq;

    #[test]
    fn test_simple_ivp_solver() {
        let df = |_t: f64, x: Vec<f64>| {vec![0.5 + 0.5*x[0].powi(2)]};

        let (ts, ys) = rk4_solver(&df, vec![3_f64.sqrt()], 0.0, 1.0, 1e-6);
        println!("{:?}, {:?}", ts[0], ts[ts.len() - 1]);
        assert_almost_eq!(ts[ts.len() - 1], 1.0);
        assert_almost_eq!(ys[ys.len() - 1][0], (0.5 + PI/3.0).tan());
    }
}
