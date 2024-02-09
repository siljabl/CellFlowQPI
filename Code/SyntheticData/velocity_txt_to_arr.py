import numpy as np
import matplotlib.pyplot as plt

input = np.loadtxt('../../Data/InitialConditions/testED.txt')

i_max = input[:,0].max()
i_min = input[:,0].min()
stepsize = input[1,1] - input[1,0]
arrsize = int((i_max - i_min) / stepsize) + 1

x_arr = np.zeros([arrsize, arrsize])
y_arr = np.zeros([arrsize, arrsize])

idx = [(input[:,0]/2).astype('int'), (input[:,1]/2).astype('int')]

x_arr[idx[0], idx[1]] = input[:,2]
y_arr[idx[0], idx[1]] = input[:,3]

np.savetxt('../../Data/InitialConditions/x_velocity.txt', x_arr)
np.savetxt('../../Data/InitialConditions/y_velocity.txt', y_arr)