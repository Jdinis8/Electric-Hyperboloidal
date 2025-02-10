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
    
    precision_id = int(f.readline())
    
    # Map precision identifier to NumPy dtype
    if precision_id == 0:
        dtype = np.float32
    elif precision_id == 1:
        dtype = np.float64
    elif precision_id == 2:
        dtype = np.longdouble
    else:
        raise ValueError(f"Unknown precision identifier: {precision_id}")
    
    num_vars, tsize = map(int, f.readline().strip().split())
    space_steps = list(map(dtype, f.readline().strip().split()))
        
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

        lowres_poisson = (list(map(dtype, f.readline().strip().split())))
        lowres_hamiltonian = (list(map(dtype, f.readline().strip().split())))
        lowres_momentum = (list(map(dtype, f.readline().strip().split())))
        lowres_z = (list(map(dtype, f.readline().strip().split())))
        
        mediumres_poisson = (list(map(dtype, f.readline().strip().split())))
        mediumres_hamiltonian = (list(map(dtype, f.readline().strip().split())))
        mediumres_momentum = (list(map(dtype, f.readline().strip().split())))
        mediumres_z = (list(map(dtype, f.readline().strip().split())))
        
        highres_poisson = (list(map(dtype, f.readline().strip().split())))
        highres_hamiltonian = (list(map(dtype, f.readline().strip().split())))
        highres_momentum = (list(map(dtype, f.readline().strip().split())))
        highres_z = (list(map(dtype, f.readline().strip().split())))
        
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
            # scatter Poisson along time
            plt.plot(space_steps, lowres_poisson, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_poisson, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_poisson, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.00001,0.00001)
            plt.xlim(0,1)
            plt.title('Poisson Constraint Convergence at time ' + str(times[j]))
            plt.savefig('python/constraintconv/poisson_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            # scatter Hamiltonian along time
            plt.plot(space_steps, lowres_hamiltonian, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_hamiltonian, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_hamiltonian, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.001,0.001)
            plt.xlim(0,1)
            plt.title('Hamiltonian Constraint Convergence at time ' + str(times[j]))
            plt.savefig('python/constraintconv/hamiltonian_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            # scatter Momentum along time
            plt.plot(space_steps, lowres_momentum, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_momentum, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_momentum, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.0005,0.0005)
            plt.xlim(0,1)
            plt.title('Momentum Constraint Convergence at time ' + str(times[j]))
            plt.savefig('python/constraintconv/momentum_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            # scatter Z along time
            plt.plot(space_steps, lowres_z, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_z, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_z, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.0001,0.0001)
            plt.xlim(0,1)
            plt.title('Z Constraint Convergence at time ' + str(times[j]))
            plt.savefig('python/constraintconv/z_conv_' + str(j) + '.png', dpi=300)
            plt.close()

with open('./output/ChargeConvergence.txt', 'r') as f:
    num_vars, tsize = map(int, f.readline().strip().split())
    space_steps = list(map(dtype, f.readline().strip().split()))
        
    #First line after this is time step
    for j in range(tsize):
        
        # Lists to store the several resolutions of the constraints
        lowres_charge = []
        mediumres_charge = []
        highres_charge = []
        
        lowres_current = []
        mediumres_current = []
        highres_current = []
        
        timestamp_line = f.readline()
        if not timestamp_line:
            break
        times.append(float(timestamp_line))

        lowres_charge = (list(map(dtype, f.readline().strip().split())))
        lowres_current = (list(map(dtype, f.readline().strip().split())))
                          
        mediumres_charge = (list(map(dtype, f.readline().strip().split())))
        mediumres_current = (list(map(dtype, f.readline().strip().split())))
        
        highres_charge = (list(map(dtype, f.readline().strip().split())))
        highres_current = (list(map(dtype, f.readline().strip().split())))
                         
        if j%1 == 0:
            # scatter Charge along time
            plt.plot(space_steps, lowres_charge, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_charge, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_charge, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.000000025,0.000000025)
            plt.xlim(0,1)
            plt.title('Charge Density at time ' + str(times[j]))
            plt.savefig('python/constraintconv/charge_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            # scatter Current along time
            plt.plot(space_steps, lowres_current, marker='o', linewidth=0.5, markersize=2, label='Low Resolution')
            plt.plot(space_steps, mediumres_current, marker='o', linewidth=0.5, markersize=2, label='Medium Resolution')
            plt.plot(space_steps, highres_current, marker='o', linewidth=0.5, markersize=2, label='High Resolution')
            plt.xlabel('Space Steps')
            plt.ylabel('Values')
            plt.legend()
            plt.ylim(-0.0000025,0.0000025)
            plt.xlim(0,1)
            plt.title('Current Density at time ' + str(times[j]))
            plt.savefig('python/constraintconv/current_conv_' + str(j) + '.png', dpi=300)
            plt.close()

if(False):            
    # Create the frames
    imgsP = glob.glob("./python/constraintconv/poisson_conv_*.png")
    imgsH = glob.glob("./python/constraintconv/hamiltonian_conv_*.png")
    imgsM = glob.glob("./python/constraintconv/momentum_conv_*.png")
    imgsZ = glob.glob("./python/constraintconv/z_conv_*.png")

    imgsCharge = glob.glob("./python/constraintconv/charge_conv_*.png")
    imgsCurrent = glob.glob("./python/constraintconv/current_conv_*.png")

    #sorting :=)
    imgsP.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
    imgsH.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
    imgsM.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
    imgsZ.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

    imgsCharge.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
    imgsCurrent.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

    ##video
    clipP = ImageSequenceClip(imgsP, fps = 20) 
    clipH = ImageSequenceClip(imgsH, fps = 20)
    clipM = ImageSequenceClip(imgsM, fps = 20)
    clipZ = ImageSequenceClip(imgsZ, fps = 20)

    clipCharge = ImageSequenceClip(imgsCharge, fps = 20)
    clipCurrent = ImageSequenceClip(imgsCurrent, fps = 20)

    clipP.write_videofile("./python/videos/video_poissonconv.mp4", fps = 30)
    clipH.write_videofile("./python/videos/video_hamiltonianconv.mp4", fps = 30)
    clipM.write_videofile("./python/videos/video_momentumconv.mp4", fps = 30)
    clipZ.write_videofile("./python/videos/video_zconv.mp4", fps = 30)
    clipCharge.write_videofile("./python/videos/video_chargeconv.mp4", fps = 30)
    clipCurrent.write_videofile("./python/videos/video_currentconv.mp4", fps = 30)