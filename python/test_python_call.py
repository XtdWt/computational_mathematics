import pytest
from computational_mathematics import (
    adaptive_quadrature,
    barycentric_lagrange_interpolation,
    bisection_method,
    chebyshev_nodes,
    cubic_spline_interpolation,
    fast_fourier_transform,
    herons_method,
    romberg_integration,
    second_derivative,
    steepest_descent,
)


def test_herons_method():
    assert herons_method(2, 2) == pytest.approx(2**0.5)


def test_bisection_method():
    temp_func = lambda x: x**2 - 2
    assert bisection_method(temp_func, 0, 2) == pytest.approx(2**0.5)


def test_bisection_method_error():
    temp_func = lambda x, y, z: x**2 - 2 + y + z
    with pytest.raises(Exception):
        bisection_method(temp_func, 0, 2)


if __name__ == "__main__":
    print(bisection_method(lambda x: x**2 - 4, -1.1, -2.1, 1000))

    p = barycentric_lagrange_interpolation([0, 2, 3], [1, 2, 4])
    for i in range(5):
        print(i, p(i))

    nodes = chebyshev_nodes(1, 4, 3)
    print(nodes)

    c = cubic_spline_interpolation([0, 1], [0, 1])
    i = 0
    while i < 1:
        print((i, c(i)))
        i += 0.1

    Xs = fast_fourier_transform([2, 3, 2, 3])
    print(Xs)

    f = lambda x: x**3

    ddf = second_derivative(f, 1, 0.001)
    print(ddf)

    print(adaptive_quadrature)

    print(romberg_integration(lambda x: x**2 + 7, -1, 4, 2))

    try:
        temp_func = lambda x, y, z: x**2 - 2 + y + z
        bisection_method(temp_func, 0, 2)
    except Exception as e:
        print(type(e))
        print(type(e.__cause__))

    print(steepest_descent)
