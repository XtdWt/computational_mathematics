use crate::rk4_solver;
use crate::root_finding::newton_raphson_method;

pub fn shooting_method<DF: Fn(f64, Vec<f64>) -> Vec<f64>>(
    df: &DF,
    ya: f64,
    yb: f64,
    t0: f64,
    tf: f64,
) -> f64 {
    let f = move |s: f64| {
        let _ts, ys = rk4_solver(
            &df,
            vec![0.0, s],
            t0,
            tf,
            1e-6,
        );
        let final_y = ys[-1];
        final_y.get(1).unwrap() - yb
    }
    newton_raphson_method(&f, &f, )
}
