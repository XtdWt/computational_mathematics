use std::cmp::Ordering;

pub trait Evaluatable {
    fn eval(&self, x: f64) -> Option<f64>;
}


pub trait Differentiable {
    fn differentiate(&self) -> Self;
}


pub trait Integrable {
    fn integrate(&self, x0: f64, y0: f64) -> Self;
}


pub fn cmp_f64(a: &f64, b: &f64) -> Ordering {
    if a < b {
        return Ordering::Less
    } else if a > b {
        return Ordering::Greater
    }
    Ordering::Equal
}
