use crate::calculus::util::t_range;


// pub fn rk4_solver(
//     df: MultivariateFunction,
//     y0: f64,
//     t0: f64,
//     tn: f64,
//     h: f64,
// ) -> (Vec<f64>, Vec<f64>) {
//     let ts = t_range(t0, tn, h);
//     let mut ys = Vec::with_capacity(ts.len());
//     ys.push(y0);
//     for i in (1..ts.len()) {
//         let y_i = ys[i-1];
//         let s1 = df(ts[i-1], y_i);
//         let s2 = df(ts[i-1] + h/2.0, ys[i-1] + (h/2.0)*s1);
//         let s3 = df(ts[i-1] + h/2.0, ys[i-1] + (h/2.0)*s2);
//         let s4 = df(ts[i-1]+h, ys[i-1] + h*s3);

//         ys.push(ys[i-1] + h * (s1 + 2.0*s2 + 2.0*s3 + s4) / 6.0)
//     }
//     return (ts, ys);
// }
