import matplotlib.pyplot as plt
import numpy as np
import numpy.linalg as linal

precision = 10e-8


def valid_input(A, f):
    A_lines     = len(A)
    A_columns   = len(A[0])

    for i in range (0, A_lines):
        if (len(A[i]) != A_columns):
            return False # A is not a matrix at all
        
    f_lines   = len(f)
    f_columns = np.array(f[0]).size

    if (A_columns != A_lines): # A is not square matrix
        return False
    
    if (float(linal.det(A)) < precision):
        return False # A - has zero determinant

    if (f_columns != 1): # f is not a vector
        return False
    
    if (f_lines != A_lines): # f and A sizes are not coordinated 
        return False
    
    return True


def build_plot(error_list, name='plot', folder='img/', dpi=500, show=False):
    
    fileName = name + '.jpg'

    plt.figure(figsize=(16/2,9/2))
    plt.plot(error_list)
    
    plt.xlim([0, len(error_list)])
    plt.ylim([0, max(error_list)])

    plt.grid(linestyle = '--', linewidth = 0.5)

    plt.xlabel("Step")
    plt.ylabel("Error value")

    plt.title(name)
    
    if show != False:
        plt.show()
    
    plt.savefig(folder + fileName, dpi=500)
    plt.clf()


def create_matrix():
    A = np.zeros((100, 100), np.float64)

    for i in range (A.shape[0]):
        A[0][i] = 1

    k = 0
    for i in range (1, A.shape[0] - 1):
        A[i][k]     = 1
        A[i][k + 1] = 10
        A[i][k + 2] = 1

        k += 1

    A[99][98] = 1
    A[99][99] = 1

    return A

def create_f():
    f = np.zeros((100), np.float64)
    for i in range (0, f.shape[0]):
        f[i] = 100 - i

    return f