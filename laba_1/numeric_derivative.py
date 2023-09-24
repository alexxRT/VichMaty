import numpy as np
import matplotlib.pyplot as plt
import pandas as pandas
import os
import math as f

result_path = "results/"

if not os.path.exists(result_path):
    os.mkdir(result_path)
    print("Result dir created successfuly!")
else:
    print("Result path dir already exist!")

#initial functions
def f_1(x):
    return f.sin(pow(x, 2))

def f_2(x):
    return f.cos(f.sin(x))

def f_3(x):
    return f.exp(f.sin(f.cos(x)))

def f_4(x):
    #range of definition for ln()
    if (x <= -3):
        return False

    return f.log(x + 3, f.e)

def f_5(x):
    #only real values are interisting
    if (x < -3): 
        return False
    
    return f.sqrt(x + 3)

f_ = [f_1, f_2, f_3, f_4, f_5]


#functions derivativies
def f_1_(x):
    return 2 * x * f.cos(pow(x, 2))

def f_2_(x):
    return -f.sin(f.sin(x)) * f.cos(x)

def f_3_(x):
    return -f.exp(f.sin(f.cos(x)))*f.cos(f.cos(x))*f.sin(x)

def f_4_(x):
    if x == -3:
        return False
    
    return 1 / (x + 3)

def f_5_(x):
    if x < -3:
        return False
    
    return 1 / (2*f.sqrt(x + 3))

f__ = [f_1_, f_2_, f_3_, f_4_, f_5_]

# h - delta parametor which is hn = (1/2ˆ(n-1)), where n = 1, 21;
# x - point of derivative
# i - index of current fuction 

def first_derivative(x, i):
    deltas = []
    for n in range (1, 22):
        h = 1 / pow(2, n - 1)

        numeric_derivativ = (f_[i](x + h) - f_[i](x)) / h
        analyt_derivative = (f__[i](x))

        deltas.append(abs(numeric_derivativ - analyt_derivative))

    return deltas

def second_derivative(x, i):
    deltas = []
    for n in range (1, 22):
        h = 1 / pow(2, n - 1)

        numeric_derivativ = (f_[i](x) - f_[i](x - h)) / h
        analyt_derivative = (f__[i](x))

        deltas.append(abs(numeric_derivativ - analyt_derivative))

    return deltas

def third_derivative(x, i):
    deltas = []
    for n in range (1, 22):
        h = 1 / pow(2, n - 1)

        numeric_derivativ = (f_[i](x + h) - f_[i](x - h)) / (2*h)
        analyt_derivative = (f__[i](x))

        deltas.append(abs(numeric_derivativ - analyt_derivative))

    return deltas

def forth_derivative(x, i):
    deltas = []
    for n in range (1, 22):
        h = 1 / pow(2, n - 1)

        numeric_derivativ = (4/3)*(f_[i](x + h)   - f_[i](x - h))   / (2*h)   \
                          - (1/3)*(f_[i](x + 2*h) - f_[i](x - 2*h)) / (4*h)   \

        analyt_derivative = (f__[i](x))

        deltas.append(abs(numeric_derivativ - analyt_derivative))

    return deltas

def fifth_derivative(x, i):
    deltas = []
    for n in range (1, 22):
        h = 1 / pow(2, n - 1)

        numeric_derivativ =  3/2  *(f_[i](x + h)   - f_[i](x - h))   / (2*h)     \
                          -  3/5  *(f_[i](x + 2*h) - f_[i](x - 2*h)) / (4*h)     \
                          +  1/10 *(f_[i](x + 3*h) - f_[i](x - 3*h)) / (6*h)     \
                          
        analyt_derivative = f__[i](x)

        delta = abs(numeric_derivativ - analyt_derivative)

        deltas.append(delta)

    return deltas

deltas_f = [first_derivative, second_derivative, third_derivative, forth_derivative, fifth_derivative]

#output to plot
deltas_  = []

# for each function set particular derivative point 
der_p = [2, 2, 2, 2, 2]

for i in range(0, 5):
    func_deltas = []
    x = der_p[i]

    for del_ in deltas_f:
        func_deltas.append(del_(x, i))

    deltas_.append(func_deltas)


def store_plots(func_num,  title):
    x_values = []
    for i in range(1, 22):
        x_values.append(1 / pow(2, i - 1))

    ax = plt.subplot()

    for i in range (0, 5):
        ax.loglog(x_values, deltas_[func_num][i])

    LegendList = ['$\\frac{f(x+h) - f(x)}{h}$',\
            
            '$\\frac{f(x) - f(x-h)}{h}$',\
            
            '$\\frac{f(x+h) - f(x-h))}{2h}$',\
            
            '$\\frac{4}{3}\\cdot\\frac{f(x+h) - f(x-h))}{2h}\
                - \\frac{1}{3}\\cdot\\frac{f(x+2h) - f(x-2h))}{4h}$',\
            
            '$\\frac{3}{2}\\cdot\\frac{f(x+h) - f(x-h))}{2h} - \\frac{3}{5}\\cdot\\frac{f(x+2h)\
                - f(x-2h))}{4h} + \\frac{1}{10}\\cdot\\frac{f(x+3h) - f(x-3h))}{6h}$'\
            ]
    
    plt.legend(LegendList, loc="upper right", fontsize="7")

    ax.set_yscale("log", base = 10)
    ax.set_xscale("log", base = 2)
    ax.set_title(title)
    ax.grid(True, "both")

    plt.ylabel("Error value")
    plt.xlabel("Step")

    plt.savefig(result_path + str(func_num) + "_derivative.png", dpi = 800, format = "png")
    plt.close()



# main
store_plots(0, "sin(xˆ2)")
store_plots(1, "cos(sin(x))")
store_plots(2, "exp(sin(cos(x)))")
store_plots(3, "ln(x + 3)")
store_plots(4, "sqrt(x + 3)")