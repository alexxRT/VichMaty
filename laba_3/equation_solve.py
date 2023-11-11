import pandas as pd
import os
import matplotlib.pyplot as plt
import time
from autograd import numpy as np, grad
from math import log
import sympy as sp

path_plot_save  = "results/equation_plot"
path_error_save = "results/error_simple/"
MAX_ITERATIONS = 100

delta_size = 0.2
precision = 0.00001
prec_num = int(np.abs(np.log10(precision)))

def PlotEquation(deltas_left, deltas_right, roots):
    x = np.arange(0.1, 10, 0.1, dtype = float)
    y = Equation(x)

    x_zero    = np.arange(0, 10, 1)
    zero_line = np.zeros(10, dtype = float)

    plt.plot(x, y)
    plt.plot(x_zero, zero_line, "--r")

    plt.title("f(x) = 2lg(x) - x/2 + 1")
    plt.xlabel("x axis")
    plt.ylabel("y axis")

    st = []

    for i in range (0, len(roots)):
        st.append("x_" + str(i + 1) + " = " + str(round(roots[i], prec_num)))

    plt.legend(st)

    y_zeros = np.zeros(len(deltas_left))
    plt.errorbar(deltas_left,  y_zeros, yerr = 0.25, fmt = ".k", ecolor="r")
    plt.errorbar(deltas_right, y_zeros, yerr = 0.25, fmt = ".k", ecolor="r")

    plt.savefig(path_plot_save, dpi = 500)

    plt.show()
    plt.clf()
    plt.close()

    return 0

def PlotErrors(res, file_name, title):
    plt.plot(res, "--r")

    plt.title(title)
    plt.xlabel("num of iteration")
    plt.ylabel("error")

    plt.savefig(path_error_save + file_name, dpi = 500)

    plt.show()
    plt.clf()
    plt.close()

def PrintResults(roots):
    print ("Results:")

    for i in range (0, len(roots)):
        output = "x_" + str(i + 1) + " = " + str(round(roots[i], prec_num))
        print(output)



def Equation(x):
    y = 2*np.log10(x) - x/2 + 1

    return y

def LogIteration(x):
    return np.abs((2*(2*np.log10(x) + 1)))

def ExpIteration(x):
    return 10**(x/4 - 1/2)

iterations = [ExpIteration, LogIteration]

def FindCompressiveIteration(x_left, x_right):
    if (x_left > x_right):
        return -1

    iter_id = -1

    for i in range (0, len(iterations)):
        iter_deriv = grad(iterations[i])

        iter_id = i
        x = x_left

        try: 
            iter_deriv(x_left)
        except(TypeError):
            # function iteration = const
            return 0.0

        while(x < x_right):
            if (np.abs(iter_deriv(x)) > 1):
                iter_id = -1
                break
            
            x += precision

        if (iter_id != -1):
            print("Iteration found! Index:", iter_id, ". For interval: [", x_left, ",", x_right, "]")
            return iterations[iter_id]

    return iter_id


# This function returns x, wich satisfies |Equation(x)| < precision
# to avoid infinite looping, when insurance = False, while makes MAX_ITERATIONS cycles
def EquationSimpleIteration(x_left, x_right, precision, insurance = False):
    #if sizes of deltas arrays are different - nothing to do
    if (len(x_left) != len(x_right)):
        return 0

    # List to store fianl results
    res = []
    x   = []

    # Iterating on initervals:
    for i in range (0, len(x_left)):
        x_i = (x_left[i] + x_right[i]) / 2

        res_i = []
        res_i.append(np.abs(Equation(x_i)))

        # depends on the interval - different iteration should be used
        iteration = FindCompressiveIteration(x_left[i], x_right[i])

        # No iteration found
        if (iteration == -1):
            print ("No appropriate iteration found")
            return 0

        if (insurance == True):
            while (res_i[len(res_i) - 1] > precision):
                x_i = iteration(x_i)

                res_i.append(np.abs(Equation(x_i)))
        else:
            count = 0
            while (count < MAX_ITERATIONS):
                x_i = iteration(x_i)

                res_i.append(np.abs(Equation(x_i)))
                count += 1

        res.append(res_i)
        x.append(x_i)
        
    return x, res

# This function returns array of localization of roots for Equation(x)
#delta_size = delta_right - delta_left
#    /\ y_axis
#    |              \     delta_right
#    |          f(x) \     ||
#    |                \    \/ 
#- - - - - - - - (- - - - -) - - - - - - - - - - - - - - -> x_axis
#    |           /\      \
#    |           ||       \ 
#            delta_left
def LocalizeZeros(x_left, x_right, delta_left, delta_right):
    if (np.abs(x_right - x_left) > delta_size):
        x_midle = (x_right + x_left) / 2

        LocalizeZeros(x_left,  x_midle, delta_left, delta_right)
        LocalizeZeros(x_midle, x_right, delta_left, delta_right)
    else: 
        if (Equation(x_left) * Equation(x_right) <= 0):
            delta_left.append(x_left)
            delta_right.append(x_right)

    return 0

def ShrinkDeltas(deltas_left, deltas_right):
    i = 0
    j = 0

    while (i < len(deltas_left)):
        j = 0
        while (j < len(deltas_right)):
            if (deltas_left[i] == deltas_right[j] and i!= j):
                new_left  = deltas_left[j]
                new_right = deltas_right[i]

                deltas_left.pop(i)
                deltas_left.pop(j)
                deltas_right.pop(i)
                deltas_right.pop(j)

                deltas_left.append (new_left)
                deltas_right.append(new_right)

                print (len(deltas_left), len(deltas_right))
            j += 1
        i += 1

def IsNutonIterationValid(x_left, x_right):
    if (x_left > x_right):
        return -1
    
    x = x_left
    equation_deriv = grad(Equation)
    second_deriv   = grad(Equation)

    while (x < x_right):
        # f'(x) != 0
        if (np.abs(equation_deriv(x)) < precision):
            return False
        
        h = 0.000000001
        limit_right = (equation_deriv(x + h) - equation_deriv(x)) / h
        limit_left  = (equation_deriv(x) - equation_deriv(x - h)) / h

        # f' does not continious
        if (np.abs(limit_right - limit_left) > precision):
            return False
        
        # f" change sign on the gieven interval
        if (second_deriv(x) * second_deriv(x + precision) < 0):
            return False
        
        x += precision

    return True

def EquationNutonIteration(x_left, x_right, precision, insurance = False):
    #if sizes of deltas arrays are different - nothing to do
    if (len(x_left) != len(x_right)):
        return 0

    # List to store fianl results
    res = []
    x   = []

    # f'(x)
    equation_deriv = grad(Equation)

    # Iterating on initervals:
    for i in range (0, len(x_left)):
        x_i = x_left[i]#(x_left[i] + x_right[i]) / 2

        res_i = []
        res_i.append(np.abs(Equation(x_i)))

        is_comressed = IsNutonIterationValid(x_left[i], x_right[i])

        if not is_comressed:
            print("Nuton iteration method is not aplicable for interval: [", x_left[i], x_right[i], "]")
            continue

        if (insurance == True):
            while (res_i[len(res_i) - 1] > precision):
                x_i = x_i - Equation(x_i) / equation_deriv(x_i)

                res_i.append(np.abs(Equation(x_i)))
        else:
            count = 0
            while (count < MAX_ITERATIONS):
                x_i = x_i - Equation(x_i) / equation_deriv(x_i)

                res_i.append(np.abs(Equation(x_i)))
                count += 1

        res.append(res_i)
        x.append(x_i)
        
    return x, res


#------------------------------------------IMPLEMENTATION---------------------------------------------

def main():

    delta_left = []
    delta_right = []

    #Now it works great!
    LocalizeZeros(0.1, 10, delta_left, delta_right)

    #Its possible for intervals to have mutural point
    #This function makes one interval from two with mutural point
    ShrinkDeltas(delta_left, delta_right)


    print ("\n\n--------------------Solving with simple iteration method--------------------\n")

    roots, res = EquationSimpleIteration(delta_left, delta_right, 0.0001, True)

    PrintResults(roots)
    PlotErrors(res[0], "method_1", "x = 10ˆ(x/4 - 1/2)")
    PlotErrors(res[1], "method_2", "x = 4lg(x) + 2")

    print("\n-----------------------------------------------------------------------------\n")

    print ("\n\n---------------------Solving with Nuton iteration method---------------------\n")

    roots, res = EquationNutonIteration(delta_left, delta_right, 0.0001, True)
    
    PrintResults(roots)
    PlotErrors(res[0], "Nuton_1", "x = 0.39")
    PlotErrors(res[1], "Nuton_2", "x = 4.68")

    print("\n-----------------------------------------------------------------------------\n")

    PlotEquation(delta_left, delta_right, roots)

if __name__ == '__main__':
    main()


