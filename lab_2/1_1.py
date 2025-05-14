import numpy as np
import matplotlib.pyplot as plt


def newton(x0, equation, equation_vp, e):
    max_iter = 1000
    counter = 0

    if (equation(x0) * equation_vp(x0) > 0):
        for _ in range(max_iter):
            counter += 1

            f_x = 3 ** x0 - 5 * x0 ** 2 + 1
            df_x = np.log(3) * 3 ** x0 - 10 * x0
            x1 = x0 - f_x / df_x

            if abs(x1 - x0) < e:
                return x1, counter

            x0 = x1

    return x0, counter

def iter(x0, e):
    max_iter = 1000
    counter = 0

    for _ in range(max_iter):
        counter += 1

        x1 = np.sqrt((3 ** x0 + 1) / 5)
        x1_pp = lambda x: (np.sqrt(5) * np.log(3) * 3 ** x) / (10 * np.sqrt(3 ** x + 1))

        q = -100
        s = -0.5

        while s < 1:
            q = max(x1_pp(s), q)
            s += 0.001

        if (q / (1 - q)) * abs(x1 - x0) <= e:
            return x1, counter

        x0 = x1

    return x0, counter

def graf(eq1, eq2):
    x = np.linspace(-1, 1, 50)
    plt.figure(figsize=(10, 6), dpi=80)

    plt.plot(x, eq1(x), label='f1(x) = 3^x')
    plt.plot(x, eq2(x), label='f2(x) = 5x^2 - 1')

    plt.legend()
    plt.title('Графики уравнений')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.show()

def main():
    e = float(input("Введите точность e: "))

    equation = lambda x: 3 ** x - 5 * x ** 2 + 1
    equation_vp = lambda x: (np.log(3) ** 2) * 3 ** x - 10
    eq1 = lambda x: 3 ** x
    eq2 = lambda x: 5 * x ** 2 - 1

    graf(eq1, eq2) 
    
    root_newton, counter_newton = newton(0.85, equation, equation_vp, e)
    root_iter, counter_iter = iter(0.85, e)

    print("\nРезультаты:")
    print(f"x по методу Ньютона: {root_newton}")
    print(f"Количество итераций: {counter_newton}")
    print(f"x по методу простых итераций: {root_iter}")
    print(f"Количество итераций: {counter_iter}")

if __name__ == "__main__":
    main()