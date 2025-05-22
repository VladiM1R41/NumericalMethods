import numpy as np
import re
import matplotlib.pyplot as plt
import math

def lagrange(f, x, test_point):
    y = [f(t) for t in x]
    assert len(x) == len(y)

    polynom_str = 'L(x) ='
    polynom_test_value = 0

    for i in range(len(x)):
        cur_enum_str = ''
        cur_enum_test = 1
        cur_denom = 1

        for j in range(len(x)):
            if i == j:
                continue
            cur_enum_str += f'(x-{x[j]:.2f})'
            cur_enum_test *= (test_point[0] - x[j])
            cur_denom *= (x[i] - x[j])

        polynom_str += f'+{(y[i] / cur_denom):.2f}' + cur_enum_str
        polynom_test_value += y[i] * cur_enum_test / cur_denom

    return polynom_str, abs(polynom_test_value - test_point[1]), polynom_test_value

def newton(f, x, test_point):
    y = [f(t) for t in x]
    assert len(x) == len(y)

    n = len(x)
    coefs = [y[i] for i in range(n)]

    for i in range(1, n):
        for j in range(n - 1, i - 1, -1):
            coefs[j] = float(coefs[j] - coefs[j - 1]) / float(x[j] - x[j - i])

    polynom_str = 'P(x) = '
    polynom_test_value = 0

    cur_multipliers_str = ''
    cur_multipliers = 1

    for i in range(n):
        polynom_test_value += cur_multipliers * coefs[i]
        
        if i == 0:
            polynom_str += f'{coefs[i]:.2f}'
        else:
            polynom_str += '+' + cur_multipliers_str + '*' + f'{coefs[i]:.2f}'

        cur_multipliers *= (test_point[0] - x[i])
        cur_multipliers_str += f'(x-{x[i]:.2f})'

    return polynom_str, abs(polynom_test_value - test_point[1]), polynom_test_value

def main():
    x_a = np.array([0, np.pi/6, 2* np.pi/6, 3 * np.pi/6])
    x_b = np.array([0, np.pi/6, np.pi/4, np.pi/2]) 
    x_star = 1.0

    equation = lambda x: np.sin(x) + x

    polynomL1, errorL1, _ = lagrange(equation, x_a, (x_star, equation(x_star)))
    polynomL2, errorL2, _ = lagrange(equation, x_b, (x_star, equation(x_star)))
    polynomN1, errorN1, _ = newton(equation, x_a, (x_star, equation(x_star)))
    polynomN2, errorN2, _ = newton(equation, x_b, (x_star, equation(x_star)))

    print("\nLagrange interpolation:")
    print(f"Points A polynom: {polynomL1}")
    print(f"Error in X*: {errorL1}")
    print(f"Points B polynom: {polynomL2}")
    print(f"Error in X*: {errorL2}\n")

    print("Newton interpolation:")
    print(f"Points A polynom: {polynomN1}")
    print(f"Error in X*: {errorN1}")
    print(f"Points B polynom: {polynomN2}")
    print(f"Error in X*: {errorN2}")

    
if __name__ == "__main__":
    main()