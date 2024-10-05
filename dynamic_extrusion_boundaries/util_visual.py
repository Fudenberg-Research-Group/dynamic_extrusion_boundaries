from scipy.ndimage.filters import gaussian_filter
import matplotlib.pylab as plt
import seaborn as sns
import numpy as np


def make_equipotential_lines(heatmap_data, line_values, num_line_x, num_line_y, label=None, legend=None,color='red', smooth_scale=1, **kwargs):

    '''
    
    heatmap_data: numpy matrix
    
    values : list of values to plot
    
    smooth_scale : how much to smooth the contour line for visualization
    
    '''
    
    d = gaussian_filter(heatmap_data.to_numpy(), sigma=smooth_scale)

    cntr1 = plt.contour(np.linspace(0, num_line_x, num_line_x * smooth_scale),

              np.linspace(0, num_line_y, num_line_y * smooth_scale),

              d, levels=line_values, colors=color, **kwargs)
    
    if label is not None:
        line, = plt.plot([], [], color=color, linestyle='dashed', alpha=0.83, label=label)

    return line