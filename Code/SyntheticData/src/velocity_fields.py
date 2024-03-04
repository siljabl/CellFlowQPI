import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp2d

# experimental velocity fields
v_file  = "testED.txt"
eps = 1e-10
   
def none(X, Y):
    u = np.zeros_like(X)
    v = np.zeros_like(Y)

    return np.array([u,v], dtype=np.float64)


def uniform(X, Y):
    u = np.ones_like(X)
    v = 0.01*np.ones_like(Y)

    return np.array([u,v], dtype=np.float64)


def curl(X, Y):
    u =  Y
    v = -X

    return np.array([u,v], dtype=np.float64)


def field_1(X, Y):
    u_exp, v_exp = field_from_file(X, Y, v_file)

    u = 2 * (0.5*v_exp + np.cos(10*Y)*(X + .3)**2 - (Y + .2)/2 + Y*np.exp(X**2))
    v = (0.5*u_exp + (X +.7) + Y) + np.sin(8*X - 2*Y)*np.sign(X)
    
    return np.array([u,v], dtype=np.float64)


def field_2(X, Y):
    u = np.sin(20*X + 3) + (10*Y-1)*np.cos(10*Y) + 10*Y
    v = np.sin(10*Y-0.4) - np.cos(10*Y-10*X) + 10*X

    return np.array([u,v], dtype=np.float64)


def normalize(u, v, u_max):
    '''
    Normalizing velocity field by max displacemet per frame
    u,v:   x,y component of velocity field
    u_max: max displacemet per frame
    '''
    scale_factor = u_max / np.sqrt(u**2 + v**2).max()

    return np.array([u, v], dtype=np.float64) * scale_factor


def field_from_file(X, Y, file):
    x, y = X[0], Y[:,0]

    up = np.loadtxt('../../Data/InitialConditions/x_velocity_' + file)
    vp = np.loadtxt('../../Data/InitialConditions/y_velocity_' + file)

    size = np.shape(up)
    xp  = np.linspace(-1, 1, size[0])
    yp  = np.linspace(-1, 1, size[1])

    u = interp2d(xp, yp, up)
    v = interp2d(xp, yp, vp)

    return np.array([u(x,y), v(x,y)], dtype=np.float64)


def cell_expansion(X, Y, loc=[0,0]):
    x = X + loc[0]
    y = Y + loc[1]
    loc = np.array(loc)
    scale = 2
    cell = norm.pdf(Y, loc=loc[1], scale=scale) * norm.pdf(X, loc=loc[0], scale=scale)
    u = cell * x
    v = cell * y

    return np.array([u,v], dtype=np.float64)