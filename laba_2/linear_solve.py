import numpy as np
import pandas as pd
import numpy.linalg as linal
import scanf as sc
from solve_methods import Gaus_method, LU_method, Yakobi_method, Zendel_method, relax_method
import io_handlers as io


#here and everywhere next:
# A - matrix of linear sysytem
# f - right part of the system
# typeA = typef = numpy array

# default value, can be specifed on start
default_precision = 5

class MODE:
    PRECISE = 1
    APPROX  = 2


def round_input(A, f, precision):
    A = np.round(A, precision)
    f = np.round(f, precision)

    return A, f


def estimate_precision(A, f, precision):
    mu = 0
    A_inverse = 0

    try:
        A_inverse = linal.inv(A)
    except linal.LinAlgError:
        print ("Matrix A does not have inverse; mu = 1")
        mu = 1

    norm_A = linal.norm(A, 2)
    norm_f = linal.norm(f, 2)
    norm_A_inverse = 0

    if (mu == 0):
        norm_A_inverse = linal.norm(A_inverse, 2)
        mu = norm_A_inverse * norm_A


    delta_A = A - np.round(A, precision)
    delta_f = f - np.round(f, precision)

    norm_delta_A = linal.norm(delta_A, 2)
    norm_delta_f = linal.norm(delta_f, 2)

    resp_delta_A = norm_delta_A / norm_A
    resp_delta_f = norm_delta_f / norm_f


    error_percent  = (mu * (1 / (1 + mu * resp_delta_A)) * (resp_delta_f + resp_delta_A)) * 100

    return error_percent

def my_func(A, f):
    return 0


def calculate_linear_straight(A, f, mode, method, input_precision = default_precision):
    if (not io.valid_input(A, f)):
        return False
    
    sol_prec = 0
    if (mode == MODE.APPROX):
        sol_prec = estimate_precision(A, f, input_precision)

        print ("Solution error estimation: ", "{:.2f}".format(sol_prec), "%")
        print ("Do u whant to improve/agree error or calculate solution with no round?\n \
               [1] - improve [2] - agree, [3] - no round")
        
        input = 0
        input = sc.scanf("%d")[0]

        while (input != 2 and input != 3):
            print("New precision:", end=" ", flush=True)
            
            input_precision = sc.scanf("%d")[0]
            sol_prec = estimate_precision(A, f, input_precision)

            print ("Solution error estimation: ", "{:.2f}".format(sol_prec), "%")

            print("[1] - improve, [2] - agree [3] - no round")
            input = sc.scanf("%d")[0]

        if (input == 3):
            return method(A, f)
        
        A_round = np.round(A, input_precision)
        f_round = np.round(f, input_precision)

        return method(A_round, f_round) 


    return method(A, f)




print ("Straight methods calculation:")

Gaus_A = io.create_matrix()
Gaus_f = io.create_f()
gaus_solution = calculate_linear_straight(Gaus_A, Gaus_f, MODE.PRECISE, Gaus_method)
print ("Gaus method error:", linal.norm(Gaus_A @ gaus_solution - Gaus_f, 2))

lu_A = io.create_matrix()
lu_f = io.create_f()
lu_solution = calculate_linear_straight(lu_A, lu_f, MODE.PRECISE, LU_method)
print ("LU method error:", linal.norm(lu_A @ lu_solution - lu_f, 2), "\n")

print("Itarative methods calculation:")
epsilon = 10e-8 #desirable precision
print("Desirable solution precision:", epsilon)

Yakobi_A = io.create_matrix()
Yakobi_f = io.create_f()
x_0 = np.zeros((Yakobi_A.shape[0]), np.float64) #starting vector
yakobi_solution, err_list_1 = Yakobi_method(io.create_matrix(), io.create_f(), x_0, epsilon)
print("Yakobi_method, Steps:", len(err_list_1))
io.build_plot(err_list_1, "Yakobi_method")

Zendel_A = io.create_matrix()
Zendel_f = io.create_f()
x_0 = np.zeros((Yakobi_A.shape[0]), np.float64) #starting vector
zendel_solution, err_list_2 = Zendel_method(Zendel_A, Zendel_f, x_0, epsilon)
print("Zendel_method, Steps:", len(err_list_2))
io.build_plot(err_list_2, "Zendel_method")

Relax_A = io.create_matrix()
Relax_f = io.create_f()
x_0 = np.zeros((Yakobi_A.shape[0]), np.float64) #starting vector
relax_solution, err_list_3 = relax_method(Relax_A, Relax_f, x_0, 0.5, epsilon)
print("Relaxation method, Steps:", len(err_list_3))
io.build_plot(err_list_3, "Relax_method")