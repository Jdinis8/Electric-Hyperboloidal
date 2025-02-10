from matplotlib import pyplot as plt

import numpy as np
from scipy.interpolate import interp1d
from scipy.interpolate import pade
from scipy.interpolate import BarycentricInterpolator

def read_data(filename):
    data = []
    with open(filename, 'rb') as file:
        # Read the precision identifier
        precision_id = np.fromfile(file, dtype=np.int32, count=1)[0]
        
        # Map precision identifier to NumPy dtype
        if precision_id == 0:
            dtype = np.float32
        elif precision_id == 1:
            dtype = np.float64
        elif precision_id == 2:
            dtype = np.longdouble
        else:
            raise ValueError(f"Unknown precision identifier: {precision_id}")
        
        # Read the number of variables
        num_variables = np.fromfile(file, dtype=np.int32, count=1)[0]
        sizex = np.fromfile(file, dtype=np.int32, count=1)[0]

        # Read the x positions
        x_positions = np.fromfile(file, dtype=dtype, count=sizex)
        
        times = []
        
        # Read the rest of the data
        while True:
            # Read the time step
            timestamp = np.fromfile(file, dtype=dtype, count=1)
            
            if timestamp.size == 0:
                break

            times.append(timestamp[0])

            # Read the variables
            variables = []
            for _ in range(num_variables):
                variable_values = np.fromfile(file, dtype=dtype, count=sizex)
                variables.append(variable_values)

            data.append((timestamp[0], variables))

    return x_positions, data, times

def extrapolate_to_x1(x, y, order=3):
    # Fit a polynomial of the given order to the data
    coeffs = np.polyfit(x, y, order)
    # Evaluate the polynomial at x=1
    return np.polyval(coeffs, 1)

def chebyshev_extrapolate(x, y, order=3):
    # Scale x to [-1,1] for better numerical stability
    x_scaled = 2 * (x - np.min(x)) / (np.max(x) - np.min(x)) - 1
    
    # Fit Chebyshev polynomials
    coeffs = np.polynomial.chebyshev.chebfit(x_scaled, y, order)
    
    # Evaluate at x=1 (which is mapped to the scaled domain)
    return np.polynomial.chebyshev.chebval(1, coeffs)

x, inputdata, times = read_data('./output/output3')

RealQNM = []
ImagQNM = []

RealQNMBefore = []
ImagQNMBefore = []

RealQNM_extrapolated = []
ImagQNM_extrapolated = []

for d in range(len(inputdata)):
    RealQNM.append(abs(inputdata[d][1][6][-1]))
    ImagQNM.append(abs(inputdata[d][1][7][-1]))
    RealQNMBefore.append(abs(inputdata[d][1][6][-2]))
    ImagQNMBefore.append(abs(inputdata[d][1][7][-2]))

for d in range(len(inputdata)):
    real_values = inputdata[d][1][6]
    imag_values = inputdata[d][1][7]
    
    #RealQNM_extrapolated.append(abs(extrapolate_to_x1(x, real_values)))
    #ImagQNM_extrapolated.append(abs(extrapolate_to_x1(x, imag_values)))
    RealQNM_extrapolated.append(abs(chebyshev_extrapolate(x, real_values)))
    ImagQNM_extrapolated.append(abs(chebyshev_extrapolate(x, imag_values)))

plt.loglog(times, RealQNM, label=r'${c_\phi}_{1-dx/2}$', color='black')
plt.loglog(times, RealQNMBefore, label=r'${c_\phi}_{1-3dx/2}$', color='black', linestyle='dashed')
plt.loglog(times, RealQNM_extrapolated, label=r'${c_\phi}_\mathscr{I}$', color='blue')
plt.ylim(1e-12,1e-4)
plt.xlim(1,times[-1])
plt.legend()
plt.title(r"Value of for $c_\phi$ at scri+")
plt.savefig('./python/Images/RealQNM.png', bbox_inches='tight')
plt.close()

plt.loglog(times, ImagQNM, label=r'${d_\phi}_{1-dx/2}$', color='red')
plt.loglog(times, ImagQNMBefore, label=r'${d_\phi}_{1-3dx/2}$', color='black', linestyle='dashed')
plt.loglog(times, ImagQNM_extrapolated, label=r'${d_\phi}_\mathscr{I}$', color='black')
plt.title(r"Value of for $d_\phi$ at scri+")
plt.legend()
plt.savefig('./python/Images/ImaginaryQNM.png', bbox_inches='tight')
plt.close()

plt.loglog(times, RealQNM, label=r'${c_\phi}_{1-dx/2}$', color='black')
plt.loglog(times, RealQNM_extrapolated, label=r'${c_\phi}_\mathscr{I}$', color='blue')
plt.loglog(times, ImagQNM, label=r'${d_\phi}_{1-dx/2}$', color='red')
plt.loglog(times, ImagQNM_extrapolated, label=r'${d_\phi}_\mathscr{I}$', color='black')
plt.legend()
plt.title(f"Comparison of $c_\phi$ and $d_\phi$ at scri+")
plt.savefig('./python/Images/RealImagQNM.png', bbox_inches='tight')
plt.close()