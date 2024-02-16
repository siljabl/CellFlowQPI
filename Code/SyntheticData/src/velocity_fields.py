import numpy as np
from scipy.stats import norm
from scipy.interpolate import interp2d

# experimental velocity fields
v_file  = "testED.txt"
eps = 1e-10
   

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


def field_1(X, Y):
    u_exp, v_exp = field_from_file(X, Y, v_file)

    u = (0.5*v_exp + np.cos(10*Y)*(X + .3)**2 - (Y + .2)/2 + Y*np.exp(X**2))
    v = (0.5*u_exp + (X +.7) + Y) + np.sin(8*X - 2*Y)*np.sign(X)
    
    return np.array([u,v], dtype=np.float64)


def field_2(X, Y):
    u = 2*(Y+1)**2 + 10*Y*np.sin(5*X)
    v = X*(X+2) + Y**3 - 6 + X + 10*Y*np.cos(X+10*Y)

    return np.array([u,v], dtype=np.float64)


def field_3(X, Y):
    u = np.sin(20*X + 3) + (10*Y-1)*np.cos(10*Y) + 10*Y
    v = np.sin(10*Y-0.4) - np.cos(10*Y-10*X) + 10*X

    return np.array([u,v], dtype=np.float64)


def field_4(X, Y):
    u_exp, v_exp = 0, 0 #field_from_file(X, Y, v_file)

    u = u_exp.T + 0.1*(2*(Y+1)**2 + 10*Y*np.sin(5*X)).T
    v = v_exp + ((X + Y*(X+0.2)) - 1)

    return np.array([u,v], dtype=np.float64)



def field_5(X,Y):
    velocity =  cell_expansion(X,Y,loc=[-.18, .2]) +\
                cell_expansion(X,Y,loc=[.1, -.28]) +\
                cell_expansion(X,Y,loc=[0.2, -.5]) +\
                cell_expansion(X,Y,loc=[.2, 0])    +\
                -cell_expansion(X,Y,loc=[-.2, -.6]) +\
                -cell_expansion(X,Y,loc=[-.5, .1]) +\
                -cell_expansion(X,Y,loc=[-.22, 0.1]) +\
                0.01*np.array([np.ones_like(X), 0*X])
    
    return velocity


def noise(X, sigma=1):
    size = np.shape(X)
    u = np.random.normal(0, sigma, size)
    v = np.random.normal(0, sigma, size)

    return np.array([u,v], dtype=np.float64)


def normalize(u, v, u_max):
    '''
    Normalizing velocity field by max displacemet per frame
    u,v:   x,y component of velocity field
    u_max: max displacemet per frame
    '''
    scale_factor = u_max / np.sqrt(u**2 + v**2).max()

    return np.array([u, v], dtype=np.float64) * scale_factor
