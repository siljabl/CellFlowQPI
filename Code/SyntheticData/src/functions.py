import numpy as np


def dI_dt(I, u):
    '''
    dI/dt according to continuity equation

    I: intensity
    u: velocity
    '''
    # Periodic boundary conditions (padding)
    pw = 1
    I_bc = np.pad(I,    mode='wrap', pad_width=pw)
    u_bc = np.pad(u[0], mode='wrap', pad_width=pw)
    v_bc = np.pad(u[1], mode='wrap', pad_width=pw)

    advection  = u_bc * np.gradient(I_bc, axis=1) + v_bc * np.gradient(I_bc, axis=0)
    divergence = I_bc * np.gradient(u_bc, axis=1) + I_bc * np.gradient(v_bc, axis=0)

    # Removing padding
    dIdt = - (advection + divergence)[pw:-pw, pw:-pw]

    return dIdt



def Euler_integration(init, u, t_steps, dt=0.1):
    '''
    Euler method

    init:    intitial conditions
    u:       velocity
    t_steps: number of iterations
    dt:      size of time step
    '''
    # Intializing empty array
    size         = np.shape(init)
    intensity    = np.zeros([t_steps + 1, *size])
    intensity[0] = init

    for t in range(t_steps):
        intensity[t+1] = intensity[t] + dI_dt(intensity[t], u) * dt

    return intensity



def RK_integration(init, u, t_steps, dt=0.1):
    '''
    Runge-Kutta method

    init:    intitial conditions
    u:       velocity
    t_steps: number of iterations
    dt:      size of time step
    '''
    # Intializing empty array
    size         = np.shape(init)
    intensity    = np.zeros([t_steps + 1, *size])
    intensity[0] = init

    for t in range(t_steps):
        k1 = dI_dt(intensity[t], u)
        k2 = dI_dt(intensity[t] + (dt/2)*k1, u)
        k3 = dI_dt(intensity[t] + (dt/2)*k2, u)
        k4 = dI_dt(intensity[t] +   dt * k3, u)

        intensity[t+1] = intensity[t] + (k1 + 2*k2 + 2*k3 + k4) * dt / 6

    return intensity



def sub_region(f, x, y):
    '''
    Generate indices of subset.
    f:   0.5 * fraction of original data that is disregarded
    x,y: arrays defining meshgrid
    '''
    # Compute indices of subset
    xmin,xmax = int( len(x) * f ), int( (len(x)-1) * (1-f) )
    ymin,ymax = int( len(y) * f ), int( (len(y)-1) * (1-f) )

    idx = [[xmin, xmax], [ymin, ymax]]

    # Define subset in meshgrid
    x_reg = [x[xmin], x[xmax], x[xmax], x[xmin], x[xmin]]
    y_reg = [y[ymin], y[xmin], y[ymax], y[ymax], y[ymin]]

    return idx, [x_reg, y_reg]