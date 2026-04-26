from computational_mathematics import (
    barycentric_lagrange_interpolation,
    bisection_method,
    chebyshev_nodes,
    cubic_spline_interpolation,
    fast_fourier_transform,
    herons_method,
    second_derivative,
)

if __name__ == "__main__":
    print(herons_method(2, 2))
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
