import numpy as np
import matplotlib.pyplot as plt

dx = np.loadtxt("output/x_displacement.txt")
dy = np.loadtxt("output/y_displacement.txt")

fig, ax = plt.subplots(1, 2, figsize=(7,3))
imx = ax[0].imshow(dx)
imy = ax[1].imshow(dy)

ax[0].set(title="x-displacement")
ax[1].set(title="y-displacement")

#fig.tight_layout()
fig.colorbar(imx, ax=ax[0])
fig.colorbar(imy, ax=ax[1])

fig.savefig("output/im_displacement.png")

I = np.loadtxt("output/intensity.txt")

fig, ax = plt.subplots(1, 1, figsize=(4,3))
im = ax.imshow(I)

fig.colorbar(im, ax=ax)
fig.savefig("output/intensity.png")

