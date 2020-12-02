def Newton(x, f, df, eps=1.0e-8, maxiter=10):
    y = x
    for n in range(maxiter):
        y = y - f(y) / df(y)
        if (abs(f(y)) < eps):
            return y
    print("収束しなかった")

def f(x):
    return x**2 - 2

def df(x):
    return 2*x

x = Newton(1.0, f, df)
print(x)
