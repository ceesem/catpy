import catmaid_interface as ci
import neuron_analytics as na
import scipy as sp
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns


def plot_neuron( nrn, **kwargs ):
    min_paths = na.minimal_paths( nrn )

    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for pth in min_paths:
        xs = [nrn.nodeloc[nd][0] for nd in pth]
        ys = [nrn.nodeloc[nd][1] for nd in pth]
        zs = [nrn.nodeloc[nd][2] for nd in pth]
        ax.plot(xs,ys,zs, **kwargs )
    ax.auto_scale_xyz([min(nrn.nodeloc.values())[0], max(nrn.nodeloc.values())[0]],
        [min(nrn.nodeloc.values())[1], max(nrn.nodeloc.values())[1]],
        [min(nrn.nodeloc.values())[2], max(nrn.nodeloc.values())[2]])
    return ax
