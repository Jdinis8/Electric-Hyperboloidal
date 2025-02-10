import numpy as np
from matplotlib import pyplot as plt
from PIL import Image
from moviepy.editor import *
import glob

def read_data(filename):
    data = []
    with open(filename, 'r') as file:
        # Read the number of variables
        num_variables = int(file.readline())

        # Read the x positions
        x_positions = [float(x) for x in file.readline().split()]

        # Read the rest of the data
        while True:
            # Read the time step
            timestamp_line = file.readline()
            if not timestamp_line:
                break
            timestamp = float(timestamp_line)

            # Read the variables
            variables = []
            for _ in range(num_variables):
                variable_line = file.readline()
                variable_values = [float(x) for x in variable_line.split()]
                variables.append(variable_values)

            data.append((timestamp, variables))
            next(file)

    return x_positions, data

x, inputdata = read_data('./output.txt')

var_index = 0

for d in range(len(inputdata)):
    if(d%10 == 0):
        plt.plot(x, [inputdata[d][1][0][i] for i in range(len(inputdata[d][1][0]))], label="Electric Field")
        plt.plot(x, [inputdata[d][1][1][i] for i in range(len(inputdata[d][1][1]))], label="Psi Field")
        plt.plot(x, [inputdata[d][1][4][i] for i in range(len(inputdata[d][1][4]))], label="Phi")
        plt.plot(x, [inputdata[d][1][5][i] for i in range(len(inputdata[d][1][5]))], label="A")
        plt.plot(x, [inputdata[d][1][6][i] for i in range(len(inputdata[d][1][6]))], label="cphi")
        plt.plot(x, [inputdata[d][1][7][i] for i in range(len(inputdata[d][1][7]))], label="dphi")
        plt.ylim(-0.015, 0.01)
        plt.xlim(0, 0.7)
        plt.legend()
        plt.title(f"Time: {inputdata[d][0]:.{7}f}")
        plt.savefig('./python/Images/func_' + str(inputdata[d][0]) + '.png')
        plt.close()

# Create the frames
frames = []
imgs = glob.glob("./python/Images/*.png")

#sorting :=)
imgs.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

##video
clip = ImageSequenceClip(imgs, fps = 20) 
clip.write_videofile("./python/output/video.mp4", fps = 30)