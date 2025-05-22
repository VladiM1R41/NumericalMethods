import numpy as np

def rectangle_trapeze(f, l, r, h, is_rectangle=True):
    if l > r:
        return None
    
    result = 0
    cur_x = l

    while cur_x < r:
        if is_rectangle:
            result += f((cur_x + cur_x + h) * 0.5)
        else:
            result += 0.5*(f(cur_x + h) + f(cur_x))
        cur_x += h

    return h*result

def simpson(f, l, r, h):
    if l > r:
        return None
    
    while ((l - r) // h) % 2 != 0:
        h *= 0.9

    result = 0
    cur_x = l + h

    while cur_x < r:
        result += f(cur_x - h) + 4 * f(cur_x) + f(cur_x + h)
        cur_x += 2 * h

    return result * h / 3

def runge_romberg(Fh, Fkh, k, p):
    return (Fh - Fkh) / (k**p - 1)

def main():
    print("Введите начальную точку x0:")
    x0 = float(input())
    print("Введите конечную точку xk:")
    xk = float(input())
    print("Введите шаг h1:")
    h1 = float(input())
    print("Введите шаг h2:")
    h2 = float(input())

    equation = lambda x: x / (x ** 3 + 8)

    rectangle_h1 = rectangle_trapeze(equation, x0, xk, h1)
    rectangle_h2 = rectangle_trapeze(equation, x0, xk, h2)
    trapeze_h1 = rectangle_trapeze(equation, x0, xk, h1, False)
    trapeze_h2 = rectangle_trapeze(equation, x0, xk, h2, False)
    simpson_h1 = simpson(equation, x0, xk, h1)
    simpson_h2 = simpson(equation, x0, xk, h2)

    rectangle_runge_rombert = runge_romberg(rectangle_h1, rectangle_h2, h2 / h1, 2)
    trapeze_runge_rombert = runge_romberg(trapeze_h1, trapeze_h2, h2 / h1, 2)
    simpson_runge_rombert = runge_romberg(simpson_h1, simpson_h2, h2 / h1, 2)

    print("\nРезультаты численного интегрирования:")
    print("\nМетод прямоугольников:")
    print(f"Шаг {h1}: {rectangle_h1:.6f}")
    print(f"Шаг {h2}: {rectangle_h2:.6f}")
    
    print("\nМетод трапеций:")
    print(f"Шаг {h1}: {trapeze_h1:.6f}")
    print(f"Шаг {h2}: {trapeze_h2:.6f}")
    
    print("\nМетод Симпсона:")
    print(f"Шаг {h1}: {simpson_h1:.6f}")
    print(f"Шаг {h2}: {simpson_h2:.6f}")
    
    print("\nОценка погрешности методом Рунге-Ромберга:")
    print(f"Для метода прямоугольников: {rectangle_runge_rombert:.6f}")
    print(f"Для метода трапеций: {trapeze_runge_rombert:.6f}")
    print(f"Для метода Симпсона: {simpson_runge_rombert:.6f}")

if __name__ == "__main__":
    main()