import numpy as np
import scipy as sc


def size(im, u):
    dim = np.shape(u)[0]
    im_size = np.shape(im)
    u_size = np.shape(u)[1:]

    size = [min(im_size[i], u_size[i]) for i in range(dim)]

    im_i = [int((im_size[i] - size[i]) / 2) for i in range(dim)]
    im_f = [size[i] + im_i[i] for i in range(dim)]

    u_i = [int((u_size[i] - size[i]) / 2) for i in range(dim)]
    u_f = [size[i] + u_i[i] for i in range(dim)]

    im_update = im[im_i[0]:im_f[0], im_i[1]:im_f[1]]
    u_update  =  u[:, u_i[0]:u_f[0], u_i[1]:u_f[1]]

    return im_update, u_update




def sub_region(f, im):
    '''
    Generate indices of subset.
    f: fraction of original data that is disregarded
    im: image data
    '''
    # Compute indices of subset
    size = np.shape(im)
    idx = [[f, size[0]-f], [f, size[1]-f]]

    # Define subset in meshgrid
    x_reg = [idx[0][0], idx[0][1], idx[0][1], idx[0][0], idx[0][0]]
    y_reg = [idx[1][0], idx[1][0], idx[1][1], idx[1][1], idx[1][0]]

    return idx, [x_reg, y_reg]




def noise(intensity):
    size = np.shape(intensity)
    norm_dist = np.random.normal(0, 1, size)
    noise = 2**8 * norm_dist / np.max(norm_dist)

    return noise.astype(int)


def smoothen(data, sigma=2):
    return sc.ndimage.gaussian_filter(data, sigma=sigma)
