import matplotlib.pyplot as plt
from moviepy.editor import *
from PIL import Image

import numpy as np
import glob
# Read the file and extract the necessary information

Maxwell = True
ScalarField = True
Einstein = True

big_error = []
small_error = []
times = []

with open('./output/pointconvergence.txt', 'r') as f:

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

    for j in range(tsize):
        big_error = []
        small_error = []
        
        times.append(float(f.readline().strip()))
        for i in range(len(space_steps)):
            big_error.append(list(map(dtype, f.readline().strip().split())))
        for i in range(len(space_steps)):
            small_error.append(list(map(dtype, f.readline().strip().split())))

        # Transpose the big_error array
        big_error_transposed = list(zip(*big_error))
        small_error_transposed = list(zip(*small_error))

        if(Maxwell): 
            # scatter the variables along time
            plt.plot(space_steps, big_error_transposed[0], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[0], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Electric Field at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/electric_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[1], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[1], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Psi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/psi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[4], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[4], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            #plt.ylim(-0.0001, 0.0001)
            plt.title(f"PC Phi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/phi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[5], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[5], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.ylim(-1e-9,1e-9)
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC A at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/A_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

        if(ScalarField):
            plt.plot(space_steps, big_error_transposed[6], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[6], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Cphi at t: {times[j]:.{7}f}")
            #plt.ylim(-0.00001, 0.00001)
            plt.savefig('python/pointconv/cphi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[7], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[7], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])            
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Dphi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/dphi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[8], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[8], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Cpi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/cpi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()

            plt.plot(space_steps, big_error_transposed[9], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[9], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Dpi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/dpi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
        
        if(Einstein):
            plt.plot(space_steps, big_error_transposed[2], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[2], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Beta at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/beta_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[3], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[3], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Alpha at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/alpha_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[10], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[10], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC trK at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/trK_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[11], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[11], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Gammarr at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/gammarr_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[12], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[12], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.ylim(-1e-11,1e-11)
            plt.title(f"PC Chi at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/chi_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[13], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[13], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Arr at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/Arr_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[14], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[14], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.ylim(-1e-8, 1e-9)
            plt.title(f"PC Lambdar at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/lambdar_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
            
            plt.plot(space_steps, big_error_transposed[15], marker='o', linewidth=0.5, markersize=2)
            plt.plot(space_steps, small_error_transposed[15], marker='o', linewidth=0.5, markersize=2)
            plt.legend(['LowRes-MedRes', '(MedRes-HighRes)*pow(f,conv)'])
            plt.xlabel('Space Steps')
            plt.ylabel('Errors')
            plt.title(f"PC Theta at t: {times[j]:.{7}f}")
            plt.savefig('python/pointconv/theta_point_conv_' + str(j) + '.png', dpi=300)
            plt.close()
if(True):
    if(Maxwell):       
        # Create the frames for electric_point_conv
        framesEl = []
        imgsEl = glob.glob("./python/pointconv/electric_point_conv_*.png")

        # Sorting the images
        imgsEl.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for electric_point_conv
        clipEl = ImageSequenceClip(imgsEl, fps=20) 
        clipEl.write_videofile("./python/videos/electric_point_conv.mp4", fps=30)

        # Create the frames for psi_point_conv
        framesPsi = []
        imgsPsi = glob.glob("./python/pointconv/psi_point_conv_*.png")

        # Sorting the images
        imgsPsi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for psi_point_conv
        clipPsi = ImageSequenceClip(imgsPsi, fps=20) 
        clipPsi.write_videofile("./python/videos/psi_point_conv.mp4", fps=30)

        # Create the frames for phi_point_conv
        framesPhi = []
        imgsPhi = glob.glob("./python/pointconv/phi_point_conv_*.png")

        # Sorting the images
        imgsPhi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for phi_point_conv
        clipPhi = ImageSequenceClip(imgsPhi, fps=20) 
        clipPhi.write_videofile("./python/videos/phi_point_conv.mp4", fps=30)

        # Create the frames for A_point_conv
        framesA = []
        imgsA = glob.glob("./python/pointconv/A_point_conv_*.png")

        # Sorting the images
        imgsA.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for A_point_conv
        clipA = ImageSequenceClip(imgsA, fps=20) 
        clipA.write_videofile("./python/videos/A_point_conv.mp4", fps=30)

    if(ScalarField):
        # Create the frames for cphi_point_conv
        framesCphi = []
        imgsCphi = glob.glob("./python/pointconv/cphi_point_conv_*.png")

        # Sorting the images
        imgsCphi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for cphi_point_conv
        clipCphi = ImageSequenceClip(imgsCphi, fps=20) 
        clipCphi.write_videofile("./python/videos/cphi_point_conv.mp4", fps=30)

        # Create the frames for dphi_point_conv
        framesDphi = []
        imgsDphi = glob.glob("./python/pointconv/dphi_point_conv_*.png")

        # Sorting the images
        imgsDphi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for dphi_point_conv
        clipDphi = ImageSequenceClip(imgsDphi, fps=20) 
        clipDphi.write_videofile("./python/videos/dphi_point_conv.mp4", fps=30)

        # Create the frames for cpi_point_conv
        framesCpi = []
        imgsCpi = glob.glob("./python/pointconv/cpi_point_conv_*.png")

        # Sorting the images
        imgsCpi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for cpi_point_conv
        clipCpi = ImageSequenceClip(imgsCpi, fps=20) 
        clipCpi.write_videofile("./python/videos/cpi_point_conv.mp4", fps=30)

        # Create the frames for dpi_point_conv
        framesDpi = []
        imgsDpi = glob.glob("./python/pointconv/dpi_point_conv_*.png")

        # Sorting the images
        imgsDpi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))

        # Create video for dpi_point_conv
        clipDpi = ImageSequenceClip(imgsDpi, fps=20) 
        clipDpi.write_videofile("./python/videos/dpi_point_conv.mp4", fps=30)
        
    if(Einstein):
        # Create the frames for alpha_point_conv
        framesBeta = []
        imgsBeta = glob.glob("./python/pointconv/beta_point_conv_*.png")
        # Sorting the images
        imgsBeta.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for alpha_point_conv
        clipBeta = ImageSequenceClip(imgsBeta, fps=20) 
        clipBeta.write_videofile("./python/videos/beta_point_conv.mp4", fps=30)
        
        # Create the frames for alpha_point_conv
        framesAlpha = []
        imgsAlpha = glob.glob("./python/pointconv/alpha_point_conv_*.png")
        # Sorting the images
        imgsAlpha.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for alpha_point_conv
        clipAlpha = ImageSequenceClip(imgsAlpha, fps=20) 
        clipAlpha.write_videofile("./python/videos/alpha_point_conv.mp4", fps=30)

        # Create the frames for trK_point_conv
        framesTrK = []
        imgsTrK = glob.glob("./python/pointconv/trK_point_conv_*.png")
        # Sorting the images
        imgsTrK.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for trK_point_conv
        clipTrK = ImageSequenceClip(imgsTrK, fps=20) 
        clipTrK.write_videofile("./python/videos/trK_point_conv.mp4", fps=30)

        # Create the frames for gammarr_point_conv
        framesGammarr = []
        imgsGammarr = glob.glob("./python/pointconv/gammarr_point_conv_*.png")
        # Sorting the images
        imgsGammarr.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for gammarr_point_conv
        clipGammarr = ImageSequenceClip(imgsGammarr, fps=20) 
        clipGammarr.write_videofile("./python/videos/gammarr_point_conv.mp4", fps=30)

        # Create the frames for chi_point_conv
        framesChi = []
        imgsChi = glob.glob("./python/pointconv/chi_point_conv_*.png")
        # Sorting the images
        imgsChi.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for chi_point_conv
        clipChi = ImageSequenceClip(imgsChi, fps=20) 
        clipChi.write_videofile("./python/videos/chi_point_conv.mp4", fps=30)

        # Create the frames for Arr_point_conv
        framesArr = []
        imgsArr = glob.glob("./python/pointconv/Arr_point_conv_*.png")
        # Sorting the images
        imgsArr.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for Arr_point_conv
        clipArr = ImageSequenceClip(imgsArr, fps=20) 
        clipArr.write_videofile("./python/videos/Arr_point_conv.mp4", fps=30)

        # Create the frames for lambdar_point_conv
        framesLambdar = []
        imgsLambdar = glob.glob("./python/pointconv/lambdar_point_conv_*.png")
        # Sorting the images
        imgsLambdar.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for lambdar_point_conv
        clipLambdar = ImageSequenceClip(imgsLambdar, fps=20) 
        clipLambdar.write_videofile("./python/videos/lambdar_point_conv.mp4", fps=30)

        # Create the frames for theta_point_conv
        framesTheta = []
        imgsTheta = glob.glob("./python/pointconv/theta_point_conv_*.png")
        # Sorting the images
        imgsTheta.sort(key=lambda entry: float(entry.rsplit('_', 1)[1].rsplit('.p',1)[0]))
        # Create video for theta_point_conv
        clipTheta = ImageSequenceClip(imgsTheta, fps=20) 
        clipTheta.write_videofile("./python/videos/theta_point_conv.mp4", fps=30)