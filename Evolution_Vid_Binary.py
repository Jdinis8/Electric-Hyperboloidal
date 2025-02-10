from matplotlib import pyplot as plt
from moviepy.editor import *
from PIL import Image

import numpy as np
import glob

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
                
        # Read the rest of the data
        while True:
            # Read the time step
            timestamp = np.fromfile(file, dtype=dtype, count=1)
            if timestamp.size == 0:
                break

            # Read the variables
            variables = []
            for _ in range(num_variables):
                variable_values = np.fromfile(file, dtype=dtype, count=sizex)
                variables.append(variable_values)

            data.append((timestamp[0], variables))

    return x_positions, data

x, inputdata = read_data('./output/output3')

if(False):
    for d in range(len(inputdata)):
        plt.plot(x, [inputdata[d][1][0][i] for i in range(len(inputdata[d][1][0]))], label=r'$E_r$', color='red')
        plt.plot(x, [inputdata[d][1][1][i] for i in range(len(inputdata[d][1][1]))], label=r'$\psi$', color='blue')
        plt.plot(x, [inputdata[d][1][2][i] for i in range(len(inputdata[d][1][2]))], label=r'$\beta_r$', color='green')
        plt.plot(x, [inputdata[d][1][3][i] for i in range(len(inputdata[d][1][3]))], label=r'$\alpha$', color='orange')
        plt.plot(x, [inputdata[d][1][4][i] for i in range(len(inputdata[d][1][4]))], label=r'$\Phi$', color='purple')
        plt.plot(x, [inputdata[d][1][5][i] for i in range(len(inputdata[d][1][5]))], label=r'$A_r$', color='black')
        plt.plot(x, [inputdata[d][1][6][i] for i in range(len(inputdata[d][1][6]))], label=r'$c_\phi$', color='pink')
        plt.plot(x, [inputdata[d][1][7][i] for i in range(len(inputdata[d][1][7]))], label=r'$d_\phi$', color='gray')
        plt.plot(x, [inputdata[d][1][8][i] for i in range(len(inputdata[d][1][8]))], label=r'$c_\pi$', color='cyan')
        plt.plot(x, [inputdata[d][1][9][i] for i in range(len(inputdata[d][1][9]))], label=r'$d_\pi$', color='magenta')
        plt.plot(x, [inputdata[d][1][10][i] for i in range(len(inputdata[d][1][10]))], label=r'$trK$', color='yellow')
        plt.plot(x, [inputdata[d][1][11][i] for i in range(len(inputdata[d][1][11]))], label=r'$\gamma_{rr}$', color='lime')
        plt.plot(x, [inputdata[d][1][12][i] for i in range(len(inputdata[d][1][12]))], label=r'$\chi$', color='teal')
        plt.plot(x, [inputdata[d][1][13][i] for i in range(len(inputdata[d][1][13]))], label=r'$A_{rr}$', color='olive')
        plt.plot(x, [inputdata[d][1][14][i] for i in range(len(inputdata[d][1][14]))], label=r'$\Lambda$', color='navy')
        plt.plot(x, [inputdata[d][1][15][i] for i in range(len(inputdata[d][1][15]))], label=r'$\Theta$', color='maroon')
        plt.ylim(-0.35, 1.1)
        #plt.ylim(-15, 15)
        #plt.ylim(-6,6)
        plt.xlim(0, 1)
        #plt.subplots_adjust(right=0.8) #This is to make room for the legend
        legend = plt.legend(
        loc='center left',
        bbox_to_anchor=(1, 0.5),
        fontsize='medium',  # Change font size
        title='Variables',  # Add a title
        title_fontsize='large',  # Title font size
        fancybox=True,  # Rounded box
        framealpha=0.5,  # Set frame transparency
        borderpad=1,  # Padding inside the legend box
        )
        # Change the legend background color
        legend.get_frame().set_facecolor('lightgrey')
        plt.title(f"Time: {inputdata[d][0]:.{7}f}")
        plt.savefig('./python/Images/func_' + str(inputdata[d][0]) + '.png', bbox_inches='tight')
        #plt.savefig(f'./python/Images/func_{ctr:03}.png', bbox_inches='tight') #for sending things to alex
        plt.close()
        
    print("Images done! Preparing video...")

    # Create the frames
    frames = []
    imgs = glob.glob("./python/Images/*.png")

    #sorting :=)
    imgs.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

    ##video
    clip = ImageSequenceClip(imgs, fps = 20) 
    clip.write_videofile("./python/videos/video.mp4", fps = 30)