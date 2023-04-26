import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


def get_dataframe_from_file(infile, header=0, index_col=0, sep='\t', max_col=21):
    '''
    Gets dataframe from asari feature table. Fixed column positions.

    Parameters
    ----------
    infile : asari feature table file
    max_col : max number of columns to read

    Returns
    -------
    Panda data frame
    '''
    return pd.read_table(infile, header=header, index_col=index_col, sep=sep,
                         usecols=range(max_col)
                         )

def asari_qc_plot(data, 
                  outfile="qc_plot.pdf",
                  height=12,
                  aspect=0.7,
                  cmap = sns.color_palette("Spectral", as_cmap=True),
                  ):
    '''
    Plot asari QC metrics in a combined figure and save to a PDF file.

    Parameters
    ----------
    data : dataframe as from asari feature table file, fixed column headers.
    outfile : output file name.
    height : height of figure, as in seaborn.
    aspect : aspect of figure, as in seaborn.
    cmap : color map, as in seaborn/matplotlib.
    
    Outputs
    -------
    A PDF file of the combined figure for QC metrics.
    '''
    sns.set_theme(style="whitegrid")
    data['log2snr'] = np.log2(data['snr'])
    data['log2area'] = np.log2(data['peak_area'])
    g = sns.relplot(
        data=data, x="log2area",
        y="goodness_fitting", 
        hue="cSelectivity",  size="log2snr",
        palette=cmap, sizes=(1, 50),
        height=height, aspect=aspect,
    )
    g.ax.xaxis.grid(True, "minor", linewidth=.25)
    g.ax.yaxis.grid(True, "minor", linewidth=.25)
    g.ax.set_ylim(0, 1.05)
    g.despine(left=True, bottom=True)
    g.savefig(outfile)
