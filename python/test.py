from computational_mathematics import (
    herons_method,
    bisection_method,
    barycentric_lagrange_interpolation,
    chebyshev_nodes,
    cubic_spline_interpolation,
)

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
