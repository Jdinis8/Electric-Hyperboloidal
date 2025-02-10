import matplotlib.pyplot as plt
from moviepy.editor import *
from PIL import Image

import glob

big_error = []
small_error = []
times = []

with open('./output/pointconvergence.txt', 'r') as f:

    num_vars, tsize = map(int, f.readline().strip().split())
    space_steps = list(map(float, f.readline().strip().split()))

    for j in range(tsize):
        big_error = []
        small_error = []
        
        times.append(float(f.readline().strip()))
        for i in range(len(space_steps)):
            big_error.append(list(map(float, f.readline().strip().split())))
        for i in range(len(space_steps)):
            small_error.append(list(map(float, f.readline().strip().split())))

        # Transpose the big_error array
        big_error_transposed = list(zip(*big_error))
        small_error_transposed = list(zip(*small_error))

        # scatter the variables along time
        plt.plot(space_steps, big_error_transposed[0], marker='o', linewidth=0.5, markersize=2)
        plt.plot(space_steps, small_error_transposed[0], marker='o', linewidth=0.5, markersize=2)
        plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
        plt.xlabel('Space Steps')
        plt.ylabel('Errors')
        plt.title(f"PC Exponential at t: {times[j]:.{7}f}")
        plt.savefig('python/pointconv/exp_point_conv_' + str(j) + '.png', dpi=300)
        plt.close()

       
# Create the frames for electric_point_conv
framesEl = []
imgsEl = glob.glob("./python/pointconv/exp_point_conv_*.png")

# Sorting the images
imgsEl.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

# Create video for electric_point_conv
clipEl = ImageSequenceClip(imgsEl, fps=20) 
clipEl.write_videofile("./python/videos/exp_point_conv.mp4", fps=30)
