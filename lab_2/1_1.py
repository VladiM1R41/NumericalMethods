import numpy as np
import matplotlib.pyplot as plt


def newton_method(x0, func, func_d, epsilon):
    max_iterations = 1000
    iteration_count = 0

    if (func(x0) * func_d(x0) > 0):
        for _ in range(max_iterations):
            iteration_count += 1

            f_x = 3 ** x0 - 5 * x0 ** 2 + 1
            df_x = np.log(3) * 3 ** x0 - 10 * x0
            x_1 = x0 - f_x / df_x

            if abs(x_1 - x0) < epsilon:
                return x_1, iteration_count

            x0 = x_1

    return x0, iteration_count

def iter_method(x0, epsilon):
    max_iterations = 1000
    iteration_count = 0

    for _ in range(max_iterations):
        iteration_count += 1

        x1 = np.sqrt((3 ** x0 + 1) / 5)
        phi_d = lambda x: (np.sqrt(5) * np.log(3) * 3 ** x) / (10 * np.sqrt(3 ** x + 1))

        q = -100
        current_x = -0.5

        while current_x < 1:
            q = max(phi_d(current_x), q)
            current_x += 0.001

        if (q / (1 - q)) * abs(x1 - x0) <= epsilon:
            return x1, iteration_count

        x0 = x1

    return x0, iteration_count

def plot_equations(func1, func2):
    x_values = np.linspace(-1, 1, 50)
    plt.figure(figsize=(10, 6), dpi=80)

    plt.plot(x_values, func1(x_values), 'b-', linewidth=2, label='$f_1(x) = 3^x$')
    plt.plot(x_values, func2(x_values), 'r-', linewidth=2, label='$f_2(x) = 5x^2 - 1$')

    plt.legend(fontsize=12, loc='upper left')
    plt.title('Графики уравнений', fontsize=14, pad=15)
    plt.xlabel('x', fontsize=12)
    plt.ylabel('y', fontsize=12)
    plt.grid(color='lightgray', linestyle=':', linewidth=0.7)
    plt.show()

def main():
    print("\nРешение нелинейного уравнения")
    print("---------------------------")
    
    epsilon = float(input("Введите точность вычислений (ε): "))

    equation = lambda x: 3 ** x - 5 * x ** 2 + 1
    equation_d = lambda x: (np.log(3) ** 2) * 3 ** x - 10
    eq1 = lambda x: 3 ** x
    eq2 = lambda x: 5 * x ** 2 - 1

    plot_equations(eq1, eq2)
    
    newton_root, newton_count = newton_method(0.85, equation, equation_d, epsilon)
    iteration_root, iteration_count = iter_method(0.85, epsilon)

    print("\nРезультаты вычислений:")
    print("---------------------")
    print(f"Метод Ньютона:")
    print(f"Корень: {newton_root:.8f}")
    print(f"Итераций: {newton_count}")
    
    print("\nМетод простых итераций:")
    print(f"Корень: {iteration_root:.8f}")
    print(f"Итераций: {iteration_count}")

if __name__ == "__main__":
    main()