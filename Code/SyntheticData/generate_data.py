import os
import sys
import argparse
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("src")
from functions import *
from velocity_fields import normalize
import integration
import plot_functions as plot
import save_functions as save

parser = argparse.ArgumentParser(description='Evolving input data according to input field')
parser.add_argument('im_data',   type=str, \
                    help='image that is used as initial condition. N.B. should be tif.')
parser.add_argument('field_dir', type=str, \
                    help='directory containing velocity field.')
parser.add_argument('-n_frames',  type=int, nargs='?', \
                    help='number of frames to generate', default=2)
parser.add_argument('-u_max',     type=int, nargs='?', \
                    help='upper limit on velocity', default=4)
parser.add_argument('-pad_width', type=int, nargs='?', \
                    help='pad width of image data',      default=60)
args = parser.parse_args()


# Input arguments
im_file   = args.im_data
field_dir = args.field_dir
n_frames  = args.n_frames
u_max     = args.u_max
pw        = args.pad_width

# Experimental parameters
d_cell = 40             # approximate cell diameter in pixels


# Creating data folder
umax_dir = field_dir + f"u_max_{str(u_max)}/"
if os.path.isdir(umax_dir) == 0:
    os.mkdir(umax_dir)

tif_dir = umax_dir + "tif/"
if os.path.isdir(tif_dir) == 0:
    os.mkdir(tif_dir)

plots_dir = umax_dir + "plots/"
if os.path.isdir(plots_dir) == 0:
    os.mkdir(plots_dir)


''' IMPORT IMAGE DATA AND VELOCITY FIELD '''
# Import initial conditions
init_cond = plt.imread(im_file)

# Pad to ensure initial conditions are well behaved
init_cond = np.pad(init_cond, mode='linear_ramp', end_values=0, pad_width=pw)


# Import velocity field and normalize to u_max
u, v = np.loadtxt(field_dir + "/x_velocity_full.txt"), np.loadtxt(field_dir + "/y_velocity_full.txt")
#if np.max(u) > u_max or np.max(v) > u_max:
u, v = normalize(*[u,v], u_max)

# Ensure intensity data and velocity field have same size.
intensity, velocity = size(init_cond, np.array([u, v]))
idx, region = sub_region(3*pw, intensity)

plot.velocity_field(velocity, intensity, region, umax_dir)
save.velocity_field(velocity, tif_dir, idx)


''' GENERATE DATA '''
# set time parameters
dt      = 0.01              # in frames
t_max   = int(n_frames)     # in frames
t_steps = int(t_max / dt)

intensity = integration.RK(intensity, velocity, t_steps=t_steps, dt=dt)
save.intensity(intensity, idx, t_max, t_steps, tif_dir, im_file)
plot.intensity(intensity, t_max, t_steps, region, umax_dir, im_file)
plot.mass_conservation(intensity, idx, dt, umax_dir)


# SAVE PARAMETERS
im_size  = np.shape(init_cond)
tif_size = (im_size[0] - 2*pw, im_size[1] - 2*pw)

config = {"im_file"         : im_file,
          "full resolution" : im_size,
          "tif resolution"  : tif_size,
          "pad_width"       : pw,
          "d_cell"          : d_cell, 
          "u_max"           : u_max, 
          "dt"              : dt}

save.config(config, umax_dir)
