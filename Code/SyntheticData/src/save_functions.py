import numpy as np
from PIL import Image

def config(config, folder):
    '''
    Saving simulation parameters
    config: dict of parameters
    folder: where config is saved
    '''
    with open(folder + "config.txt", 'w') as f: 
        f.write('Units\nlenght: pixels\ntime: frames\n\n') 

        for key, value in config.items():  
            f.write('%s:%s\n' % (key, value))



def velocity_field(velocity, idx, folder):
    '''
    Saving velocity field as txt
    velocity: x and y components of velocity field
    idx:      indices of subset
    folder:   where data is saved
    '''
    xlim, ylim = idx
    u = velocity[0, xlim[0]:xlim[1], ylim[0]:ylim[1]]
    v = velocity[1, xlim[0]:xlim[1], ylim[0]:ylim[1]]
   
    np.savetxt(folder + "x_velocity.txt", u)
    np.savetxt(folder + "y_velocity.txt", v)



def intensity(intensity, idx, t_max, t_steps, folder, file):
    '''
    Saving intensities as tif
    intensity: intensities [t, s]
    idx:       indices of subset
    t_max:     number of frames
    t_steps:   iterations
    folder:    where data is saved
    file:      filename
    '''
    xlim, ylim = idx

    for frame in range(t_max+1):
        # frame to index
        i = int(frame * t_steps / t_max)
        I = intensity[i, xlim[0]:xlim[1], ylim[0]:ylim[1]]

        # Save frame as tif
        I_tif  = I.astype(np.uint8)
        im_tif = Image.fromarray(I_tif)
        im_tif.save(folder + "tif/" + file + f"_%00i.tif" % (frame))