#!/usr/bin/python3

"""
Combriat T. 2024
This small script input original images and DIC results and output the original images with the DIC 
results superimposed on top.
These results are plotted as points, which position is updated according to the DIC fields.
This is useful to visualize the performances of the DIC.

"""




#################
# library imports
#################
import matplotlib.pyplot as plt  # plotting
import numpy as np  # numerics
import math as math  # math
import skvideo.io  # video opener
import skimage.io as io
import skimage.io
from tqdm import tqdm  # fancy progress bar
from scipy.optimize import curve_fit  # fitting
import os as os # folder and file manipulation


input_folder = "RescaledRegistered_sm"
prename = "Well1-1_resc_reg"
postname = ".tif"

try:
    os.mkdir(input_folder+os.sep+"DIC_super")
except:
    pass

N_min = 2
N_max = 100

N_pts = 50


img = io.imread(input_folder+os.sep+prename+str(N_min)+postname)

pts_X, pts_Y = np.meshgrid(np.linspace(0,len(img[0]),N_pts,endpoint=False),np.linspace(0,len(img),N_pts,endpoint=False))



plt.figure(figsize=(15,15))
for i in tqdm(range(N_min,N_max+1)):
    img = io.imread(input_folder+os.sep+prename+str(i)+postname)
    if (i==N_min):
        d_X = np.zeros_like(img)
        d_Y = np.zeros_like(img)
    else:
        d_X = io.imread("delta_X_av_6"+os.sep+str(i)+".tif")
        d_Y = io.imread("delta_Y_av_6"+os.sep+str(i)+".tif")

        #print(d_X.shape)
        #print(pts_X.shape)

        pts_X, pts_Y = np.meshgrid(np.linspace(0,len(img[0]),N_pts,endpoint=False),np.linspace(0,len(img),N_pts,endpoint=False))
        pts_X = pts_X+d_X[np.array(pts_Y,dtype=int),np.array(pts_X,dtype=int)]
        to_rm = np.logical_or(pts_X > len(img[0])-2,pts_X<0)
        to_rm2 = np.logical_or(pts_Y > len(img)-2,pts_Y<0)
        to_rm = np.logical_or(to_rm,to_rm2)
        pts_X[to_rm] = 0
        pts_Y[to_rm] = 0
        pts_Y = pts_Y+d_Y[np.array(pts_Y,dtype=int),np.array(pts_X,dtype=int)]

        to_rm = np.logical_or(pts_X > len(img[0])-2,pts_X<0)
        to_rm2 = np.logical_or(pts_Y > len(img)-2,pts_Y<0)
        to_rm = np.logical_or(to_rm,to_rm2)
        pts_X[to_rm] = 0
        pts_Y[to_rm] = 0

        """
        for i in range(len(pts_X)):
            print(i)
            pts_X[i] += d_X[np.array(pts_Y[i],dtype=int),np.array(pts_X[i],dtype=int)]
        """
    plt.xlim((0,len(img[0])))
    plt.ylim((0,len(img)))             
    plt.imshow(img)
    plt.scatter(pts_X,pts_Y,c="r",alpha=0.5,s=3)
    #plt.show()
    plt.savefig(input_folder+os.sep+"DIC_super"+os.sep+str(i)+".png")
    plt.clf()
    
