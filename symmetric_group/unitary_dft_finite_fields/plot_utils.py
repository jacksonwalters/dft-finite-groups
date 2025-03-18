"""
plot the discrete log of a matrix valued in a finite field
since log(0) = -infinity, we set this to -1, and color it black

- `F`: is the field
- `M`: is a matrix of discrete log values of elements of F
- `path`: the prefix of the directory to save the file
- `title`: the title of the plot, used in the saved name replacing ` ` with `_`
- `normalize`: optional parameter to rescale values by F.order() - 1, or take log1p of values

"""
def plot_discrete_log(F, M, path, title, normalize=None):
    import numpy as np
    import matplotlib.pyplot as plt
    from matplotlib.colors import ListedColormap, BoundaryNorm

    if normalize == "rescale":
        rescaled_M = M / (F.order() - 1)  # Normalize values to [0, 1]
        value_range = 1
    if normalize == "log1p":
        M_numeric = np.array(M, dtype=np.float64) #convert array to numpy floats
        log1p_M = np.log1p(M_numeric) #apply log1p transformation
        pos_log1p_M = log1p_M[(log1p_M != -np.inf) & (log1p_M != 0)] #mask out the -np.inf values and zero values for normalization
    if normalize == None:
        value_range = F.order() - 1
    
    cmap = plt.cm.hsv
    num_colors = min(256, 256)
    new_colors = np.vstack(([0, 0, 0, 1], cmap(np.linspace(0, 1, num_colors))))
    custom_cmap = ListedColormap(new_colors)

    # Normalize: exclude -np.inf values for the min/max scaling
    if normalize == "log1p":
        norm = BoundaryNorm([np.min(pos_log1p_M), *np.linspace(np.min(pos_log1p_M), np.max(pos_log1p_M), num_colors)], custom_cmap.N)
    else:
        norm = BoundaryNorm([-1] + list(np.linspace(0, value_range, num_colors)), custom_cmap.N)

    if normalize == "rescale":
        plt.imshow(rescaled_M, cmap=custom_cmap, norm=norm, interpolation="nearest")
    if normalize == "log1p":
        plt.imshow(log1p_M, cmap=custom_cmap, norm=norm, interpolation="nearest")
    if normalize == None:
        plt.imshow(M, cmap=custom_cmap, norm=norm, interpolation="nearest")
    
    plt.title(title, fontsize=16)
    plt.colorbar()
    
    # Adjust filename based on normalization
    plot_path_filename = f"{path}{title.replace(' ', '_')}_norm={normalize}.png"
    
    plt.savefig(plot_path_filename, dpi=300, bbox_inches="tight")
    plt.show()

#plot the complexified version of the uDFT matrix over a finite field
from sage.plot.matrix_plot import matrix_plot
from sage.functions.other import arg
def plot_arg_complex(U_complex, title):
    U_arg = U_complex.apply_map(lambda x: arg(x))  # find the argument of each element
    plot = matrix_plot(U_arg, cmap='hsv', colorbar=True, title=title)  # plot the matrix
    filename = "plots/dft_matrix/" + title.replace(" ", "_") + ".png"
    plot.save(filename, dpi=300)  # Save the plot as a PNG file with high resolution
    return plot