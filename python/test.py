from comp_math import herons_method, bisection_method

print(herons_method(2, 2))
print(bisection_method(lambda x: x**2 - 4, -1.1, -2.1, 1000))