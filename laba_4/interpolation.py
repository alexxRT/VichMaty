import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

year_data = [92228496,
             106021537,
             123202624,
             132164569,
             151325798,
             179323175,
             203211926,
             226545805,
             248709873,
             281421906]

true_val = 308745538

ref_data = np.array(year_data, dtype=float)

# P(x) = f(x_0) + f(x_0, x_1)(x - x_0) + f(x_0, x_1, x_2)(x - x_0)(x - x_1) ....

def NutonPolynom(devided_diferences, x_data, x):
    if len(devided_diferences) == 0:
        print("Not able to evaluate P(x), because lack of coefficients!")
        print("Devided diferences size:", len(devided_diferences))

        return 0

    # Value of Nuton Plynom in a given point x
    P_x = 0

    for i in range (0, len(devided_diferences)):
        brackets_mul = 1
        for j in range (0, i):
            brackets_mul *= (x - x_data[j])

        brackets_mul *= devided_diferences[i]
        P_x += brackets_mul

    return P_x

# data - values to be interpolated
# x    - argumants where data(x[i]) = data[i]
def NutonInterpol(data: list, x: list):
    # coefficients in P(x) described above
    f = np.zeros(len(data))

    devided_differences = []

    diff_0 = data.copy()
    devided_differences.append(diff_0)

    for i in range (0, len(data) - 1):
        diff = []
        for j in range (0, len(devided_differences[i]) - 1):
            dev_diff = (devided_differences[i][j + 1] - devided_differences[i][j])             \
                                        / (x[j  + i + 1] - x[j])
            
            diff.append(dev_diff)

        devided_differences.append(diff)
    
    for i in range (0, len(devided_differences)):
        f[i] = devided_differences[i][0]
    

    return f

def DisplayNutonInterpol(data, x):
    f = NutonInterpol(data, x.tolist())

    plt.scatter(x, data)

    x_interpol = np.arange(min(x), max(x), 0.1)
    y_interpol = np.zeros(len(x_interpol))

    for i in range(0, len(x_interpol)):
        y_interpol[i] = NutonPolynom(f, x, x_interpol[i])

    plt.plot(x_interpol, y_interpol, "r--")
    plt.xlabel("Year")
    plt.ylabel("Population quantity in USA")
    plt.title("Intorpolated data with Nuton Polynom")
    plt.show()
    plt.close()

    predicted_val = NutonPolynom(f, x, 2010)
    print("\nNuton Polynom predicted:", predicted_val)
    print("Method error", round((np.abs(predicted_val - true_val) / true_val) * 100, 2), "%\n")


def f_fitted(x, a, b):
    return a*np.exp(x*b)

def DisplayMLSInterpol(data, x):
    data = data / 10e7
    x = x / 10e2

    params, pcov = curve_fit(f_fitted, x, data)

    x_interpol = np.arange(x.min(), x.max(), 0.01, dtype=float)
    y_interpol = np.zeros(len(x_interpol), dtype=float)

    for i in range(0, len(x_interpol)):
        y_interpol[i] = f_fitted(x_interpol[i], *params)

    plt.scatter(x, data)
    plt.plot(x_interpol, y_interpol, "r--")
    plt.xlabel("Year")
    plt.ylabel("Population quantity in USA")
    plt.title("Intorpolated data with MLS, y = a*exp(b*x)")
    plt.show()
    plt.close()

    predicted_val = f_fitted(2.01, *params)
    true_scaled = true_val / 10e7
    print("\nMLS spline predicted:", predicted_val)
    print("Method error", round((np.abs(predicted_val - true_scaled) / true_scaled) * 100, 2), "%\n")

def EvaluateA_k(data, x):
    a = np.zeros(len(x) - 1)
    a = data[0:len(a)]

    return a

def EvaluateC_k(data, x):
    # boundary conditions c_0, c_n = 0
    c = np.zeros(len(x))

    # we have n+1 points -> n intervals -> n-1 unknown c values
    A_size = len(x)
    A = np.zeros((A_size, A_size), dtype=float)

    f = np.zeros(A_size, dtype=float)

    h = np.zeros(len(x) - 1)
    for i in range (0, len(x) - 1):
        h[i] = x[i + 1] - x[i]

    for i in range (0, len(A[0]) - 2):
        A[i][i]     = h[i]
        A[i][i + 1] = 2*(h[i] + h[i + 1])
        A[i][i + 2] = h[i + 1]

        f[i] = 3*(((data[i + 2] - data[i + 1]) / h[i + 1]) - ((data[i + 1] - data[i]) / h[i]))

    # Now we need to drop first and last columns + last two lines
    A = np.delete(A, A_size - 1, 0)
    A = np.delete(A, A_size - 2, 0)

    A = np.delete(A, 0, 1)
    A = np.delete(A, A_size - 2, 1)

    f = np.delete(f, A_size - 1, 0)
    f = np.delete(f, A_size - 2, 0)

    c = np.linalg.solve(A, f)
    c = np.insert(c, 0, 0.)
    c = np.append(c, 0.)

    return c

def EvaluateB_k(data, x, c):
    b = np.zeros(len(x) - 1)
    for i in range (0, len(b)):
        b[i] = ((data[i + 1] - data[i]) / (x[i + 1] - x[i])) - \
        (x[i + 1] - x[i])*(2*c[i] + c[i + 1]) / 3

    return b

def EvaluateD_k(data, x, c):
    d = np.zeros(len(x) - 1)
    for i in range (0, len(d)):
        d[i] = (c[i + 1] - c[i]) / (3*(x[i + 1] - x[i]))

    return d

def QubicPolynom(x, x_k, a, b, c, d):
    return a + b*(x - x_k) + c*(x - x_k)**2 + d*(x - x_k)**3


def DisplayQubicSplineInterpol(data, x):
    a = EvaluateA_k(data, x)
    c = EvaluateC_k(data, x)
    b = EvaluateB_k(data, x, c)
    d = EvaluateD_k(data, x, c)

    plt.scatter(x, data)

    for i in range (0, len(a)):
        x_interpol = np.arange(x[i], x[i + 1], 0.1)
        y_interpol = np.zeros(len(x_interpol))

        for j in range (0, len(x_interpol)):
            y_interpol[j] = QubicPolynom(x_interpol[j], x[i], a[i], b[i], c[i], d[i])

        plt.plot(x_interpol, y_interpol, "r--")

    plt.xlabel("Year")
    plt.ylabel("Population quantity in USA")
    plt.title("Intorpolated data with classic Qubic Spline")

    plt.show()
    plt.close()

    predicted_val = QubicPolynom(2010, x[len(x) - 2], a[len(a) - 1], b[len(b) - 1], c[len(c) - 2], d[len(d) - 1])
    print("\nQubic spline predicted:", predicted_val)
    print("Method error", round((np.abs(predicted_val - true_val) / true_val) * 100, 2), "%\n")


if __name__ == "__main__":
    x_n = np.arange(1910, 2010, 10, dtype=float)

    DisplayNutonInterpol(ref_data.copy(), x_n.copy())
    DisplayMLSInterpol(ref_data.copy(), x_n.copy())
    DisplayQubicSplineInterpol(ref_data.copy(), x_n.copy())