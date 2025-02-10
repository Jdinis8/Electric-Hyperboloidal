import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
from moviepy.editor import *
import glob

times = []

def read_data(filename):
    data = []
    data2 = []
    data3 = []
    data4 = []
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
        sizex = np.fromfile(file, dtype=np.int32, count=1)[0]

        # Read the x positions
        x_positions = np.fromfile(file, dtype=dtype, count=sizex)

        # Read the rest of the data
        while True:
            # Read the time step
            timestamp = np.fromfile(file, dtype=dtype, count=1)
            if timestamp.size == 0:
                break
            times.append(timestamp[0])
            
            # Read the variables
            poisson_values = np.fromfile(file, dtype=dtype, count=sizex)
            Hr_values = np.fromfile(file, dtype=dtype, count=sizex)
            Mr_values = np.fromfile(file, dtype=dtype, count=sizex)
            Zr_values = np.fromfile(file, dtype=dtype, count=sizex)

            data.append((timestamp[0], poisson_values))
            data2.append((timestamp[0], Hr_values))
            data3.append((timestamp[0], Mr_values))
            data4.append((timestamp[0], Zr_values))

    return x_positions, data, data2, data3, data4

x, poisson, hr, mr, zr = read_data('./output/Constraints')
      
for d in range(len(poisson)):
    plt.plot(x, [poisson[d][1][i] for i in range(len(poisson[d][1]))], label="Poisson Constraint")
    plt.plot(x, [hr[d][1][i] for i in range(len(hr[d][1]))], label="Hamiltonian Constraint")
    plt.plot(x, [mr[d][1][i] for i in range(len(mr[d][1]))], label="Momentum Constraint")
    plt.plot(x, [zr[d][1][i] for i in range(len(zr[d][1]))], label="Z Constraint")
    plt.ylim(-0.00000025, 0.00000025)
    #plt.ylim(-0.3, 0.3)
    plt.xlim(0, 1)
    plt.legend()
    plt.title("Constraints at time t = " + str(times[d]))
    plt.savefig('./python/Images/constraints_' + str(poisson[d][0]) + '.png')
    plt.close()

print("Images done! Preparing video...")

# Create the frames
frames = []
imgs = glob.glob("./python/Images/constraints_*.png")

#sorting :=)
imgs.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

##video
clip = ImageSequenceClip(imgs, fps = 20) 
clip.write_videofile("./python/videos/constraints.mp4", fps = 30)