import matplotlib.pyplot as plt
import numpy as np 
from PIL import Image
from moviepy.editor import *
import glob

# Read the file and extract the necessary information

big_error = []
small_error = []
times = []

with open('./output/pointconvergence.txt', 'r') as f:
    num_vars, tsize = map(int, f.readline().strip().split())
    space_steps = list(map(float, f.readline().strip().split()))
    
    #First line after this is time step
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
        #Plotting evertyhing summed

        sum_big_error_transposed = np.sum(big_error_transposed, axis=0)
        sum_small_error_transposed = np.sum(small_error_transposed, axis=0)

        if j%5 == 0:
            # scatter the variables along time
            plt.plot(space_steps, sum_big_error_transposed, marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, sum_small_error_transposed, marker='o', linewidth=0.5, markersize=2)
            plt.legend(['Sum of Big Errors', 'Sum of Small Errors'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            #plt.ylim(-0.0001,0.0001)
            #plt.xlim(0,1)
            plt.title(f"Point Wise Sum of Errors at time {times[j]:.{5}f}")
            plt.savefig('python/pointconv/sum_errors_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
# Create the frames
frames = []
imgs = glob.glob("./python/pointconv/*.png")

#sorting :=)
imgs.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

##video
clip = ImageSequenceClip(imgs, fps = 20) 
clip.write_videofile("./python/videos/video_pointconv.mp4", fps = 30)
