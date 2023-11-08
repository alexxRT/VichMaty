import matplotlib.pyplot as plt
from autograd import numpy as np


def BuildSystemPlot(x_vals, y_vals, file_vals, file_sys):
    plt.xlabel("Iteration")
    plt.ylabel("(x, y) value")
    plt.title("cos(x - 1) + y = 0.5\nx - cosy = 3")
    plt.xticks(np.arange(1, len(x_vals) + 1), fontsize = 6)

    plt.plot(np.arange(1, len(x_vals) + 1), x_vals, marker='.', markersize=4)
    plt.plot(np.arange(1, len(y_vals) + 1), y_vals, marker='.', markersize=4)

    LegendList = "x", "y"
    plt.legend(LegendList, loc="upper right", fontsize="10")
    plt.savefig(file_vals, dpi=500)
    
    plt.show()
    plt.close()

    plt.title("System\ncos(x - 1) + y = 0.5\nx - cosy = 3")
    plt.xlabel("x axis")
    plt.ylabel("y axis")

    x = np.arange(-10, 10, 0.1, dtype = float)
    y = 0.5 - np.cos(x - 1)
    plt.plot(x, y)

    y = np.arange(-10, 10, 0.1, dtype = float)
    x = 3 + np.cos(y)
    plt.plot(x, y)

    plt.savefig(file_sys, dpi=500)

    plt.show()
    plt.close

    plt.clf()


def ReverseJacobian(x, y):
    f1x = -np.sin(x - 1)
    f1y = -1 
    f2x = 1
    f2y = np.sin(y)

    J = np.array([[f1x, f1y], [f2x, f2y]])

    rev_J = np.linalg.inv(J)
    return rev_J


def System(x, y):
    f1 = np.cos(x - 1) - y - 0.5
    f2 = x - np.cos(y) - 3
    return [f1, f2]

def SolveNonLinearSystem(precision):
    x = [2, 1]
    x_vals = []
    y_vals = []

    while np.abs(System(x[0], x[1])[0] - System(x[0], x[1])[1]) > precision:
        x = x - ReverseJacobian(x[0], x[1]) @ System(x[0], x[1])
        x_vals.append(x[0])
        y_vals.append(x[1])

    
    BuildSystemPlot(x_vals, y_vals, "results/system_solution_plot.png", "results/system_plot.png")
    print("System solution = ", x)


def main():
    precision = 0.00001
    SolveNonLinearSystem(precision)

if __name__ == '__main__':
    main()