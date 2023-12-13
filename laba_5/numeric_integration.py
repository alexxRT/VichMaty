import numpy as np

def IntegralFunction(x):
    if x == 0:
        return 1.

    return np.sin(x) / x

def IntegrateTrapezoid(f, x, h):
    I_val = 0

    for i in range (0, len(x)):
        cur  = h*i
        next = h*(i + 1)

        if (next < len(x)):
            step = x[next] - x[cur]
            s_h = step*(f[next] + f[cur])/2

            I_val += s_h

    return I_val

def IntegrateSimpson(f, x):
    I_val = 0

    for i in range (0, len(x), 2):
        cur    = i
        middle = i + 1
        next   = i + 2

        if (next < len(x)):
            step = x[next] - x[cur]
            s_h = step*(f[next] + 4*f[middle] + f[cur])/6
            I_val += s_h

    return I_val

def RichardsonExtropol(I_2h, I_h, p):
    return I_h + (I_h - I_2h) / (2**p - 1)


x = [0, 0.25, 0.5, 0.75, 1.0, 1.25, 1.5, 1.75, 2]
x_data = np.array(x,  dtype = float)
f_data = np.array([], dtype = float)

for x_i in x_data:
    f_data = np.append(f_data, IntegralFunction(x_i))

I_2h = IntegrateTrapezoid(f_data, x_data, 2)
I_h  = IntegrateTrapezoid(f_data, x_data, 1)
I_S  = IntegrateSimpson(f_data, x_data) 

print(f"I_h  = {I_h}")
print(f"I_2h = {I_2h}")
print(f"I_S  = {I_S}")
print(f"I_R  = {RichardsonExtropol(I_2h, I_h, 2)}")