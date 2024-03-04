import numpy as np
import matplotlib.pyplot as plt


def im_velocity(velocity, dir):
    u_norm = np.sqrt(velocity[0]**2 + velocity[1]**2)
    u_min  = np.min(u_norm)
    u_max  = np.max(u_norm)

    title = [r'$u_x$', r'$u_y$']

    fig, ax = plt.subplots(1,2, figsize=(8,3.5), sharey=True)
    for i in range(2):
        ax[i].set(title=title[i])
                  #xlabel=r'$x ~[pixels]$', \
                  #ylabel=r'$y ~[pixels]$')
        
        im = ax[i].imshow(velocity[i])#, vmin=u_min, vmax=u_max)
        fig.colorbar(im, ax=ax[i]) # not sure this works
        
    fig.tight_layout()
    fig.savefig(dir + "xy_velocity.png")



def velocity_field(velocity, init_cond, region, dir):
    size  = np.shape(init_cond)
    I_min = np.min(init_cond)
    I_max = np.max(init_cond)

    u_norm = np.sqrt(velocity[0]**2 + velocity[1]**2)
    u_min  = np.min(u_norm)
    u_max  = np.max(u_norm)
    u_flip = -np.array([np.flip(velocity[i], axis=0) * (-1)**i for i in range(2)])

    x   = np.arange(size[0])
    y   = np.arange(size[1])
    X,Y = np.meshgrid(x,y)

    fig, ax = plt.subplots(1,1, figsize=(6,5))
    ax.set(xlabel=r'$x ~[pixels]$', \
           ylabel=r'$y ~[pixels]$', \
           title=r'Velocity field, $|\bar{u}(x,y)|\in$' + f' [{u_min:0.0f}, {u_max:0.0f}]')
    
    ax.streamplot(X, Y, *u_flip, density=2, linewidth=1, color='firebrick')
    ax.plot(region[0],  region[1], alpha=1,   ls="dashed", color="k", lw=3)
    im=ax.imshow(init_cond, vmin=I_min, vmax=I_max, cmap="Greys_r", alpha=0.65)
    
    fig.colorbar(im)
    fig.savefig(dir + "plots/velocity_field.png")



def intensity(intensity, t_max, t_steps, region, dir, im_file):
    I_min = np.min(intensity[0])
    I_max = np.max(intensity[0])

    filename = im_file.split('/')[-1]
    filename = filename.split('.')[0]

    for frame in range(t_max+1):
        i = int(frame * t_steps / t_max)

        # Plot and save as png
        fig, ax = plt.subplots(1, 1)
        ax.set(title=f"Frame: {frame:0.0f}",)# \
               #xlabel=r'$x ~/~ d_{cell}$', \
               #ylabel=r'$y ~/~ d_{cell}$')
        im=ax.imshow(intensity[i], 
                     origin="upper", vmin=I_min, vmax=I_max)


        ax.plot(region[0], region[1], lw=3, color="k", alpha=1, ls="dashed")
        
        fig.colorbar(im)
        fig.savefig(dir + "plots/" + filename + f"_%00i.png" % (frame), format = "png")



def mass_conservation(intensity, idx, dt, dir):
    xlim, ylim = idx

    total_mean = [np.mean(I) for I in intensity]
    tif_mean   = [np.mean(I[xlim[0]:xlim[1], ylim[0]:ylim[1]]) for I in intensity]
    time = [i * dt for i in range(len(total_mean))]

    fig, ax = plt.subplots(1,1)
    ax.plot(time, total_mean  / total_mean[0],  '.', label=r"$\langle I \rangle_{total}$")
    ax.plot(time, tif_mean / tif_mean[0], '.', label=r"$\langle I \rangle_{tif}$")
    ax.set(xlabel="Frame", ylabel="I(t) / I(0)")
    ax.legend()
    
    fig.tight_layout()
    fig.savefig(dir + "plots/mean_intensity.png")

    return fig