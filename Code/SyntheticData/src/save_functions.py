import numpy as np
from PIL import Image

def config(config, dir):
    '''
    Saving simulation parameters
    config: dict of parameters
    dir: where config is saved
    '''
    with open(dir + "/config.txt", 'w') as f: 
        f.write('Units\nlenght:pixels\ntime:frames\n\n') 

        for key, value in config.items():  
            f.write('%s:%s\n' % (key, value))



def velocity_field(velocity, dir, idx=[[0,-1], [0,-1]]):
    '''
    Saving velocity field as txt
    velocity: x and y components of velocity field
    idx:      indices of subset
    dir:   where data is saved
    '''
    xlim, ylim = idx
    u = velocity[0, xlim[0]:xlim[1], ylim[0]:ylim[1]]
    v = velocity[1, xlim[0]:xlim[1], ylim[0]:ylim[1]]
    
    if idx == [[0,-1], [0,-1]]:
        np.savetxt(dir + "/x_velocity_full.txt", u)
        np.savetxt(dir + "/y_velocity_full.txt", v)

    else:
        np.savetxt(dir + "/x_velocity_tif.txt", u)
        np.savetxt(dir + "/y_velocity_tif.txt", v)



def intensity(intensity, idx, t_max, t_steps, tif_dir, im_file):
    '''
    Saving intensities as tif
    intensity: intensities [t, s]
    idx:       indices of subset
    t_max:     number of frames
    t_steps:   iterations
    dir:       where data is saved
    file:      filename
    '''
    xlim, ylim = idx
    filename = im_file.split('/')[-1]
    filename = filename.split('.')[0]

    for frame in range(t_max+1):
        # frame to index
        i = int(frame * t_steps / t_max)
        I = intensity[i, xlim[0]:xlim[1], ylim[0]:ylim[1]]

        # Save frame as tif
        I[I < 0] = 0    # avoid negative overflows
        I_tif  = I.astype(np.uint8)
        im_tif = Image.fromarray(I_tif)
        im_tif.save(tif_dir + filename + f"_%00i.tif" % (frame))