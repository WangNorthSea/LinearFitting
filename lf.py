#encoding:UTF-8
import numpy as np
import matplotlib.pyplot as plt

# 输入数据
x = np.array([2, 3, 4, 5, 6])
y = np.array([4.725, 7.1, 9.475, 11.85, 14.275])

# 最小二乘法求拟合直线斜率、截距
def Least_squares(x, y):
    x_ = x.mean()
    y_ = y.mean()
    m = np.zeros(1)
    n = np.zeros(1)
    k = np.zeros(1)
    p = np.zeros(1)
    for i in range(0, x.size - 1):
        k = (x[i]-x_)* (y[i]-y_)
        m += k
        p = np.square( x[i]-x_ )
        n = n + p
    a = m / n
    b = y_ - a * x_
  
    return a, b

# 求可决系数R方
def square_R(x, y, a, b):
    # Caution: f = ax + b
    x_ = x.mean()
    y_ = y.mean()
    square_x_ = (x * x).mean()
    square_y_ = (y * y).mean()
    xy_ = (x * y).mean()

    up = xy_ - (x_ * y_)
    down = ((square_x_ - (x_ * x_)) * (square_y_ - (y_ * y_))) ** 0.5

    square_R0 = (up / down) * (up / down)
    square_R = round(square_R0, 8)
    return square_R

# 求y、x、b、a标准差
def standardDeviation(x, y, a, b):
    # Caution: f = bx + a!!!
    n = x.size
    f = b * x + a
    x_ = x.mean()
    square_x_ = (x * x).mean()
    
    square_sig_y_up = 0
    for i in range(0, y.size - 1):
        square_sig_y_up += (y[i] - f[i]) ** 2
    square_sig_y_down = n - 2
    sigma_y0 = (square_sig_y_up / square_sig_y_down) ** 0.5

    sigma_a0 = ((square_x_ / (n * (square_x_ - (x_ * x_)))) ** 0.5) * sigma_y0
    
    sigma_b0 = ((1 / n * (square_x_ - (x_ * x_))) ** 0.5) * sigma_y0

    square_sig_x_up = 0
    for i in range(0, x.size - 1):
        square_sig_x_up += (x[i] - (f[i] / b - a / b)) ** 2
    square_sig_x_down = n - 2
    sigma_x0 = (square_sig_x_up / square_sig_x_down) ** 0.5

    sigma_x = round(sigma_x0, 8)
    sigma_y = round(sigma_y0, 8)
    sigma_a = round(sigma_a0, 8)
    sigma_b = round(sigma_b0, 8)

    return sigma_x, sigma_y, sigma_a, sigma_b

if __name__ == '__main__':
    a, b = Least_squares(x, y)
    square_R = square_R(x, y, a, b)
    sigma_x, sigma_y, sigma_b, sigma_a = standardDeviation(x, y, b, a)
    print a, b, square_R, sigma_x, sigma_y
    y1 = a * x + b
    plt.figure(figsize=(10, 5), facecolor='w')
    plt.plot(x, y, 'ro', lw=1, markersize=6)
    plt.plot(x, y1, 'r-', lw=1, markersize=6)
    plt.grid(b=True, ls=':')
    plt.xlabel(u'I(mA)', fontsize=16)
    plt.ylabel(u'U(mV)', fontsize=16)
    plt.show()
