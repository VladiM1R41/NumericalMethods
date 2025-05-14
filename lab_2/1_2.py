import numpy as np
import matplotlib.pyplot as plt


def system_equations(x, a):
    return np.array([
        x[0] - np.cos(x[1]) - a,
        x[1] - np.sin(x[0]) - a
    ])

def jacobian_matrix(x):
    return np.array([
        [1, np.sin(x[1])],
        [-np.cos(x[0]), 1]
    ])

def newton_method(x0, a, epsilon):
    iteration_count = 0

    for _ in range(1000):
        iteration_count += 1

        delta = np.linalg.solve(jacobian_matrix(x0), 
                      -system_equations(x0, a))
        x0 += delta

        if np.linalg.norm(delta) < epsilon:
            break

    return x0, iteration_count

def iter_method(x0, a, epsilon):
    iteration_count = 0
    current_x1, current_x2 = x0
    
    while True:
        iteration_count += 1
        new_x1 = a + np.cos(current_x2)
        new_x2 = a + np.sin(current_x1)

        q = 0.5

        if (q / (1 - q)) * \
           np.sqrt((new_x1 - current_x1) ** 2 + (new_x2 - current_x2) ** 2) <= epsilon:
            return np.array([new_x1, new_x2]), iteration_count
        
        current_x1, current_x2 = new_x1, new_x2

def plot_equations(a, plot_range):
    x_values = np.linspace(*plot_range, 400)
    
    x1_curve = np.cos(x_values) + a
    x2_curve = np.sin(x_values) + a

    plt.figure(figsize=(10, 8))
    
    plt.plot(x1_curve, x_values, 'b-', linewidth=2, label='$x_1 = \cos(x_2) + a$')
    plt.plot(x_values, x2_curve, 'r-', linewidth=2, label='$x_2 = \sin(x_1) + a$')

    plt.xlabel('Ось $x_1$', fontsize=12)
    plt.ylabel('Ось $x_2$', fontsize=12)
    plt.axhline(0, color='black', linewidth=0.8)
    plt.axvline(0, color='black', linewidth=0.8)
    plt.grid(color='lightgray', linestyle=':', linewidth=0.7)
    plt.legend(fontsize=10, loc='upper right')
    plt.title('Визуализация системы уравнений\n$x_1 - \cos(x_2) = a$ и $x_2 - \sin(x_1) = a$', 
              fontsize=14, pad=20)
    plt.show()

def main():
    print("\nРешение системы нелинейных уравнений")
    print("----------------------------------")
    
    epsilon = float(input("Введите точность вычислений (ε): "))
    a = float(input("Введите параметр системы (a): "))

    plot_equations(a, plot_range=(-4, 4))

    x0 = np.array([2., 2.])

    newton_solution, newton_count = newton_method(x0, a, epsilon)
    iteration_solution, iteration_count = iter_method(x0, a, epsilon)

    print("\nРезультаты вычислений:")
    print("---------------------")
    print(f"Метод Ньютона:")
    print(f"x₁ = {newton_solution[0]}")
    print(f"x₂ = {newton_solution[1]}")
    print(f"Количество итераций: {newton_count}")
    
    print("\nМетод простых итераций:")
    print(f"x₁ = {iteration_solution[0]}")
    print(f"x₂ = {iteration_solution[1]}")
    print(f"Количество итераций: {iteration_count}")

if __name__ == "__main__":
    main()