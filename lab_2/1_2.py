import numpy as np
import matplotlib.pyplot as plt


def F(x, a):
    return np.array([
        x[0] - np.cos(x[1]) - a,
        x[1] - np.sin(x[0]) - a
    ])

def J(x):
    return np.array([
        [1, np.sin(x[1])],
        [-np.cos(x[0]), 1]
    ])

def newton(x0, a, e):
    counter = 0

    for _ in range(1000):
        counter += 1

        delta = np.linalg.solve(J(x0), -F(x0, a))
        x0 += delta

        if np.linalg.norm(delta) < e:
            break

    return x0, counter

def iter(x0, a, epsilon):
    counter = 0
    x1, x2 = x0
    
    while True:
        counter += 1
        x1_new = a + np.cos(x2)
        x2_new = a + np.sin(x1)

        q = 0.5

        if (q / (1 - q)) * np.sqrt((x1_new - x1) ** 2 + (x2_new - x2) ** 2) <= epsilon:
            return np.array([x1_new, x2_new]), counter
        
        x1, x2 = x1_new, x2_new

def graf(a, x_range):
    x = np.linspace(*x_range, 400)
    
    x1_curve1 = np.cos(x) + a
    x2_curve2 = np.sin(x) + a

    plt.figure(figsize=(8, 6))
    
    plt.plot(x1_curve1, x, label='$x_1 = \cos(x_2) + a$')
    plt.plot(x, x2_curve2, label='$x_2 = \sin(x_1) + a$')


    plt.xlabel('$x_1$')
    plt.ylabel('$x_2$')
    plt.axhline(0, color='black', linewidth=0.5)
    plt.axvline(0, color='black', linewidth=0.5)
    plt.grid(color='gray', linestyle='--', linewidth=0.5)
    plt.legend()
    plt.title('Graphs of the system: $x_1 - \cos(x_2) = a$ and $x_2 - \sin(x_1) = a$')
    plt.show()

def main():
    e = float(input("Введите точность e: "))
    a = float(input("Введите параметр a: "))

    graf(a, x_range=(-4, 4))

    x0_1 = np.array([2., 2.])

    root_newton, counter_newton = newton(x0_1, a, e)
    root_iter, counter_iter = iter(x0_1, a, e)

    print("\nРезультаты:")
    print(f"x1 по методу Ньютона: {root_newton[0]}")
    print(f"x2 по методу Ньютона: {root_newton[1]}")
    print(f"Количество итераций (Ньютон): {counter_newton}")
    print(f"\nx1 по методу простых итераций: {root_iter[0]}")
    print(f"x2 по методу простых итераций: {root_iter[1]}")
    print(f"Количество итераций (итерации): {counter_iter}")

if __name__ == "__main__":
    main()