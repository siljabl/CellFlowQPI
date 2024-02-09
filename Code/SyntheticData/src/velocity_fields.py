import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp2d

# infinitesimal value to avoid overflow
eps    = 1e-8
   

def cell_expansion(X, Y, loc=[0,0]):
    x = X + loc[0]
    y = Y + loc[1]
    loc = np.array(loc)
    scale = 2
    cell = norm.pdf(Y, loc=loc[1], scale=scale) * norm.pdf(X, loc=loc[0], scale=scale)
    u = cell * x
    v = cell * y

    return np.array([u,v])


def field_1(X, Y):
    u = X**2 - 2*(Y+1)**2 - (X-1)
    v = (X + Y) - 1
    
    return np.array([u,v])


def field_2(X, Y):
    u = 2*(Y+1)**2 + 10*Y*np.sin(X)
    v = X*(X+2) + Y**3 - 6 + X + 10*Y*np.cos(X)

    return np.array([u,v])


def field_3(X, Y):
    u = np.sin(2*X + 3) + (Y-1)*np.cos(Y) + Y
    v = np.sin(Y-0.4) - np.cos(Y-X) + X

    return np.array([u,v])


def field_4(X,Y):
    velocity =  cell_expansion(X,Y,loc=[-1.8, 2]) +\
                cell_expansion(X,Y,loc=[1, -2.8]) +\
                cell_expansion(X,Y,loc=[0.2, -5]) +\
                cell_expansion(X,Y,loc=[2, 0])    +\
                -cell_expansion(X,Y,loc=[-2, -6]) +\
                -cell_expansion(X,Y,loc=[-.5, 1]) +\
                -cell_expansion(X,Y,loc=[-2.2, 0.1]) +\
                0.01*np.array([np.ones_like(X), 0*X])
    
    return velocity


def field_exp(x, y, d_cell):
    up = np.loadtxt('../../Data/InitialConditions/x_velocity.txt')
    vp = np.loadtxt('../../Data/InitialConditions/y_velocity.txt').T

    size = np.shape(up)

    xp  = len(x) * np.linspace(-.5, .5, size[0]) / d_cell
    yp  = len(y) * np.linspace(-.5, .5, size[1]) / d_cell

    u = interp2d(xp, yp, up)
    v = interp2d(xp, yp, vp)

    return np.array([u(x,y), v(x,y)])


def noise(X, sigma=1):
    size = np.shape(X)
    u = np.random.normal(0, sigma, size)
    v = np.random.normal(0, sigma, size)

    return np.array([u,v])


def normalize(u, v, u_max):
    '''
    Normalizing velocity field by max displacemet per frame
    u,v:   x,y component of velocity field
    u_max: max displacemet per frame
    '''
    scale_factor = u_max / np.sqrt(u**2 + v**2).max()

    return np.array([u, v]) * scale_factor
