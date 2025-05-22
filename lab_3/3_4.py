import numpy as np
import matplotlib.pyplot as plt


def df(x, y, x_):
    assert len(x) == len(y)

    for interval in range(len(x)):
        if x[interval] <= x_ < x[interval+1]:
            i = interval
            break

    a1 = (y[i+1] - y[i]) / (x[i+1] - x[i])
    a2 = ((y[i+2] - y[i+1]) / (x[i+2] - x[i+1]) - a1) / (x[i+2] - x[i]) * (2*x_ - x[i] - x[i+1])

    return a1 + a2

def d2f(x, y, x_):
    assert len(x) == len(y)

    for interval in range(len(x)):
        if x[interval] <= x_ < x[interval+1]:
            i = interval
            break

    num = (y[i+2] - y[i+1]) / (x[i+2] - x[i+1]) - (y[i+1] - y[i]) / (x[i+1] - x[i])

    return 2 * num / (x[i+2] - x[i])

def main():
    

    x = np.array([-1.0, -0.4, 0.2, 0.6, 1.0])
    y = np.array([-1.4142, -0.55838, 0.27870, 0.84008, 1.4142])
    
    x_ = 0.2

    first_derivative = df(x, y, x_)
    second_derivative = d2f(x, y, x_)

    print("\nРезультаты вычислений:")
    print(f"Первая производная f'({x_}) = {first_derivative:.4f}")
    print(f"Вторая производная f''({x_}) = {second_derivative:.4f}")


if __name__ == "__main__":
    main()