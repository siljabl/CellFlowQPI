import numpy as np


def size(im, u):
    dim = np.shape(u)[0]
    im_size = np.shape(im)
    u_size = np.shape(u)[1:]

    size = [min(im_size[i], u_size[i]) for i in range(dim)]

    im_i = [int((im_size[i] - size[i]) / 2) for i in range(dim)]
    im_f = [size[i] - im_i[i] for i in range(dim)]

    u_i = [int((u_size[i] - size[i]) / 2) for i in range(dim)]
    u_f = [size[i] - u_i[i] for i in range(dim)]

    im_update = im[im_i[0]:im_f[0], im_i[1]:im_f[1]]
    u_update  =  u[:, 0:size[0], 0:size[1]]

    return im_update, u_update




def sub_region(f, im):
    '''
    Generate indices of subset.
    f: fraction of original data that is disregarded
    im: image data
    '''
    # Compute indices of subset
    size = np.shape(im)
    xmin,xmax = int( size[0] * f ), int( (size[0]-1) * (1-f) )
    ymin,ymax = int( size[1] * f ), int( (size[1]-1) * (1-f) )

    idx = [[xmin, xmax], [ymin, ymax]]

    # Define subset in meshgrid
    x_reg = [xmin, xmax, xmax, xmin, xmin]
    y_reg = [ymin, xmin, ymax, ymax, ymin]

    return idx, [x_reg, y_reg]