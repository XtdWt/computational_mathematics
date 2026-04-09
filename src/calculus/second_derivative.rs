use crate::Function;
use crate::calculus::util::{DerivativeType};


fn second_order_forward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (2.0*function(x) - 5.0*function(x+h) + 4.0*function(x+2.0*h) - function(x+3.0*h)) / (h.powf(2.0))
}


fn second_order_backward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (2.0*function(x) - 5.0*function(x-h) + 4.0*function(x-2.0*h) - function(x-3.0*h)) / (h.powf(2.0))
}


fn second_order_central_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (function(x + h) - 2.0 * function(x) + function(x - h)) / (h.powf(2.0))
}


pub fn second_derivative(
    function: Function,
    x: f64,
    h: f64,
    method: DerivativeType,
) -> f64 {
    match method {
        DerivativeType::Forward => second_order_forward_difference(function, x, h),
        DerivativeType::Backward => second_order_backward_difference(function, x, h),
        DerivativeType::Central => second_order_central_difference(function, x, h),
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn second_derivative_forward() {}

    #[test]
    fn second_derivative_backward() {}

    #[test]
    fn second_derivative_central() {}
}