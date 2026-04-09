use crate::Function;


pub enum DerivativeType {
    Forward,
    Backward,
    Central
}


fn first_order_forward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (-3.0 * function(x) + 4.0 * function(x+h) - function(x+2.0*h)) / (2.0 * h)
}


fn first_order_backward_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (3.0 * function(x) - 4.0 * function(x-h) + function(x-2.0*h)) / (2.0 * h)
}


fn first_order_central_difference(
    function: Function,
    x: f64,
    h: f64,
) -> f64 {
    (function(x + h) - function(x - h)) / (2.0 * h)
}


pub fn first_derivative(
    function: Function,
    x: f64,
    h: f64,
    method: DerivativeType,
) -> f64 {
    match method {
        DerivativeType::Forward => first_order_forward_difference(function, x, h),
        DerivativeType::Backward => first_order_backward_difference(function, x, h),
        DerivativeType::Central => first_order_central_difference(function, x, h),
    }
}


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
    fn first_derivative_forward() {

    }

    #[test]
    fn first_derivative_backward() {

    }

    #[test]
    fn first_derivative_central() {

    }

    #[test]
    fn second_derivative_forward() {

    }

    #[test]
    fn second_derivative_backward() {

    }

    #[test]
    fn second_derivative_central() {

    }
}