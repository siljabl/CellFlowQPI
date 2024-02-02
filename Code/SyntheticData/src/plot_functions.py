import numpy as np
import matplotlib.pyplot as plt


def velocity_field(velocity, X, Y, init_cond, region, folder):
    I_min = np.min(init_cond)
    I_max = np.max(init_cond)

    u_norm = np.sqrt(velocity[0]**2 + velocity[1]**2)
    u_min = np.min(u_norm)
    u_max = np.max(u_norm)

    fig, ax = plt.subplots(1,1)
    ax.set(xlim=(X[0,0], X[0,-1]),
           ylim=(X[-1,0], X[-1,-1]),
           xlabel=r'$x ~/~ d_{cell}$', 
           ylabel=r'$y ~/~ d_{cell}$', 
           title=r'Velocity field, $|\bar{u}(x,y)|\in$' + f' [{u_min:0.0f}, {u_max:0.0f}]')
    
    #ss = 8
    #ax.quiver(X[0:-1:ss, 0:-1:ss], Y[0:-1:ss, 0:-1:ss], velocity[0, 0:-1:ss, 0:-1:ss], velocity[1,0:-1:ss, 0:-1:ss], width=0.002, color='firebrick')
    ax.streamplot(X, Y, *velocity, density=4, linewidth=1, color='firebrick')
    ax.plot(region[0], region[1], lw=3, color="k", alpha=1, ls="dashed")
    im=ax.imshow(init_cond, origin="lower", extent=[X[0,0], X[0,-1], X[-1,0], X[-1,-1]], 
                 vmin=I_min, vmax=I_max, cmap="Greys_r", alpha=0.65)
    
    fig.colorbar(im)
    fig.savefig(folder + "velocity_field.pdf")



def intensity(intensity, X, t_max, t_steps, region, folder, file):
    I_min = np.min(intensity[0])
    I_max = np.max(intensity[0])

    for frame in range(t_max+1):
        i = int(frame * t_steps / t_max)

        # Plot and save as png
        fig, ax = plt.subplots(1, 1)
        ax.set(title=f"Frame: {frame:0.0f}", 
               xlabel=r'$x ~/~ d_{cell}$',
               ylabel=r'$y ~/~ d_{cell}$')
        im=ax.imshow(intensity[i], extent=[X[0,0], X[0,-1], X[-1,0], X[-1,-1]], 
                     origin="upper", vmin=I_min, vmax=I_max)


        ax.plot(region[0], region[1], lw=3, color="k", alpha=1, ls="dashed")
        
        fig.colorbar(im)
        fig.savefig(folder + file + f"_%00i.png" % (frame), format = "png")



def mass_conservation(intensity, idx, dt, folder):
    xlim, ylim = idx

    total_mean  = [np.mean(I) for I in intensity]
    region_mean = [np.mean(I[xlim[0]:xlim[1], ylim[0]:ylim[1]]) for I in intensity]
    time = [i * dt for i in range(len(total_mean))]

    fig, ax = plt.subplots(1,1)
    ax.plot(time, total_mean  / total_mean[0],  '.', label=r"$\langle I \rangle_{total}$")
    ax.plot(time, region_mean / region_mean[0], '.', label=r"$\langle I \rangle_{subset}$")
    ax.set(xlabel="Frame", ylabel="I(t) / I(0)")
    ax.legend()
    
    fig.tight_layout()
    fig.savefig(folder + "mean_intensity.png")

    return fig