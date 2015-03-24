'''
!python bam_to_qualdepth.py subsample.bam  | cat > f
run graph_qualdepth f -o subsampled.png
'''
import sys
sys.path.append('/home/AMED/michael.panciera/projects/ngs_mapper/ngs_mapper/')
import graph_qualdepth

import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec


def simple_plot_depths(depths_array):

    xvals = range(len(depths_array))
    yvals = depths_array
    maxdepth = max(depths_array) + 100
    color='green'
    title='simple_plot'

    fig = plt.figure()
    fig.set_size_inches( 20.0, 8.0 )
    gs = gridspec.GridSpec(1,2, width_ratios=[20,1])
    ax = plt.subplot(gs[0]) 
    graph_qualdepth.plot_depths(ax, xvals, depths_array,
            maxdepth, color, title)

    outputfile = 'foo.png'
    fig.savefig( outputfile, bbox_inches='tight', dpi=100, pad_inches=0.1 )
    fig.show()


if __name__ == '__main__':
    pass



