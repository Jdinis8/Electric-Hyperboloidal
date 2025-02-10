import matplotlib.pyplot as plt
import numpy as np 
from PIL import Image
from moviepy.editor import *
import glob

# Read the file and extract the necessary information

times = []

# fh - the amount of times we increase the step size
# p - the order of convergence of our code

fh = 3
p = 4

with open('./output/ConstraintConvergence.txt', 'r') as f:
    num_vars, tsize = map(int, f.readline().strip().split())
    space_steps = list(map(float, f.readline().strip().split()))
        
    #First line after this is time step
    for j in range(tsize):
        
        # Lists to store the several resolutions of the constraints
        lowres_poisson = []
        lowres_hamiltonian = []
        lowres_momentum = []
        lowres_z = []

        mediumres_poisson = []
        mediumres_hamiltonian = []
        mediumres_momentum = []
        mediumres_z = []

        highres_poisson = []
        highres_hamiltonian = []
        highres_momentum = []
        highres_z = []
        
        timestamp_line = f.readline()
        if not timestamp_line:
            break
        times.append(float(timestamp_line))

        lowres_poisson = (list(map(float, f.readline().strip().split())))
        lowres_hamiltonian = (list(map(float, f.readline().strip().split())))
        lowres_momentum = (list(map(float, f.readline().strip().split())))
        lowres_z = (list(map(float, f.readline().strip().split())))
        
        mediumres_poisson = (list(map(float, f.readline().strip().split())))
        mediumres_hamiltonian = (list(map(float, f.readline().strip().split())))
        mediumres_momentum = (list(map(float, f.readline().strip().split())))
        mediumres_z = (list(map(float, f.readline().strip().split())))
        
        highres_poisson = (list(map(float, f.readline().strip().split())))
        highres_hamiltonian = (list(map(float, f.readline().strip().split())))
        highres_momentum = (list(map(float, f.readline().strip().split())))
        highres_z = (list(map(float, f.readline().strip().split())))
        
        # Multiply the elements of mediumres arrays by pow(f,p)
        mediumres_poisson = [x * pow(fh, p) for x in mediumres_poisson]
        mediumres_hamiltonian = [x * pow(fh, p) for x in mediumres_hamiltonian]
        mediumres_momentum = [x * pow(fh, p) for x in mediumres_momentum]
        mediumres_z = [x * pow(fh, p) for x in mediumres_z]

        # Multiply the elements of highres arrays by pow(f, 2*p)
        highres_poisson = [x * pow(fh, 2*p) for x in highres_poisson]
        highres_hamiltonian = [x * pow(fh, 2*p) for x in highres_hamiltonian]
        highres_momentum = [x * pow(fh, 2*p) for x in highres_momentum]
        highres_z = [x * pow(fh, 2*p) for x in highres_z]
                                   
        if j%5 == 0:
            # scatter Hamiltonian along time
            plt.plot(space_steps, lowres_hamiltonian, linestyle=':', color='black', label='Low Resolution')
            plt.plot(space_steps, mediumres_hamiltonian, linestyle='-.', color='red', label='Medium Resolution')
            plt.plot(space_steps, highres_hamiltonian, color='#315c31', label='High Resolution')
            plt.xlabel('r')
            plt.ylabel('H')
            plt.legend()
            plt.ylim(-0.007,0.007)
            plt.xlim(0,1)
            plt.title(f'Hamiltonian Constraint Convergence at time {times[j]:{5}f}')
            plt.savefig('python/constraintconv/hamiltonian_conv_' + str(j) + '.pdf')
            plt.close()
            
            # scatter Momentum along time
            plt.plot(space_steps, lowres_momentum, linestyle=':', color='black', label='Low Resolution')
            plt.plot(space_steps, mediumres_momentum, linestyle='-.', color='red', label='Medium Resolution')
            plt.plot(space_steps, highres_momentum, color='#315c31', label='High Resolution')
            plt.xlabel('r')
            plt.ylabel('M')
            plt.legend()
            plt.ylim(-0.005,0.005)
            plt.xlim(0,1)
            plt.title(f'Momentum Constraint Convergence at time {times[j]:{5}f}')
            plt.savefig('python/constraintconv/momentum_conv_' + str(j) + '.pdf')
            plt.close()
            
            # scatter Z along time
            plt.plot(space_steps, lowres_z, linestyle=':', color='black', label='Low Resolution')
            plt.plot(space_steps, mediumres_z, linestyle='-.', color='red', label='Medium Resolution')
            plt.plot(space_steps, highres_z, color='#315c31', label='High Resolution')
            plt.xlabel('r')
            plt.ylabel('Z')
            plt.legend()
            plt.ylim(-0.00005,0.00005)
            plt.xlim(0,1)
            plt.title(f'Z Constraint Convergence at time {times[j]:{5}f}')
            plt.savefig('python/constraintconv/z_conv_' + str(j) + '.pdf')
            plt.close()