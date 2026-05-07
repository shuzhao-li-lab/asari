import numpy as np
import matplotlib.pyplot as plt
from asari.tools.ms2 import extract_all_spectra_form_file

def summarize_ms_file(infile, nspec_plot=50, offset_n=100, width=10, height=6, 
                      title='', outfile='ms_summary.pdf'):
    '''Plot characteristics and summary of a mass spec file, GC/LC-MS or MS/MS.
    {
        number_spectra_level1
        number_spectra_level2
        max_rtime
    }
    '''
    [ms1_spectra, ms2_spectra, others] = extract_all_spectra_form_file(
        infile, min_intensity=100, MS2_peak_limit=999999)
    # the spectra extraction is sequential by scan number, e.g. 'sp1'
    plot_msn(ms1_spectra, ms2_spectra, 
             nspec_plot=nspec_plot, offset_n=offset_n, 
             width=width, height=height, title=title, outfile=outfile)
    
    
def plot_msn(ms1_spectra, ms2_spectra, 
             nspec_plot=50, offset_n=100, width=10, height=6, title='', outfile=None):
    
    def get_median_height(peaks):  # use median instead
        return np.log10(np.median([x[1] for x in peaks]))
    def get_max_height(peaks):
        z = [x[1] for x in peaks]
        if z: return max(z)
        else: return 0
    
    number_spectra_level1, number_spectra_level2 = len(ms1_spectra), len(ms2_spectra)
    max_rtime = max([x['rtime'] for x in ms1_spectra+ms2_spectra])
    if number_spectra_level1 + number_spectra_level2 < nspec_plot + offset_n:
        NN = range(nspec_plot)
    else:
        NN = range(offset_n, offset_n+nspec_plot)
    # select scans
    to_use_ids = ['sp' + str(ii+1) for ii in NN]
    ms1_plotdata, ms2_plotdata = [], []
    for x in ms1_spectra:
        if x['id'] in to_use_ids:
            x['level'] = 1
            x['num_peaks'] = len(x['peaks'])
            x['median_intensity'] = get_median_height(x['peaks'])
            ms1_plotdata.append(x)
    for x in ms2_spectra:
        if x['id'] in to_use_ids:
            x['level'] = 2
            x['num_peaks'] = len(x['peaks'])
            x['median_intensity'] = get_median_height(x['peaks'])
            ms2_plotdata.append(x)
    # sort by rtime; no need
    # to_use.sort(key=lambda x: x['rtime'])
    text_offset = 0.1
    _rotation_ = 0
    if not ms2_spectra:     # then MS1 labels get crowded
        _rotation_ = 60
            
    fig, ax = plt.subplots(2, 1, figsize=(width, height))
    fig.suptitle(title, fontsize=12, x=0.05)
    
    # top panel for scan illustration
    xx = [x['rtime'] for x in ms1_plotdata]
    yy = [x['median_intensity'] for x in ms1_plotdata]
    ax[0].vlines(xx,
               [2] * len(ms1_plotdata), # ymin
               yy,   # ymax
               colors='b', linestyles='solid', label='MS1'
               )
    for x in ms1_plotdata:
        ax[0].text(x['rtime'], x['median_intensity'] + text_offset, 
                   x['num_peaks'], 
                   rotation=_rotation_,
                   fontsize=6, ha='center')

    ax[0].vlines([x['rtime'] for x in ms2_plotdata],
               [2] * len(ms2_plotdata), # ymin
               [x['median_intensity'] for x in ms2_plotdata],   # ymax
               colors='r', linestyles='dashed', 
               label='MS2'
               )
    for x in ms2_plotdata:
        ax[0].text(x['rtime'], x['median_intensity'] + text_offset, x['num_peaks'], fontsize=6, ha='center')
        
    ax[0].set_ylim(ymax=max(yy) + 0.5)
    ax[0].set_xlabel("Retention time (sec)")
    ax[0].set_ylabel("Median log10 intensity")
    ax[0].legend(loc='upper right', ncols=2, bbox_to_anchor=(1, 1.2))
    ax[0].spines['top'].set_visible(False)
    ax[0].spines['right'].set_visible(False)
    ax[0].set_title("Example scans", fontsize=12)
    
    # bottom panel for XIC
    ax[1].plot(
        [x['rtime'] for x in ms1_spectra], 
        [get_max_height(x['peaks']) for x in ms1_spectra], 'b-'
    )
    note = '\n'.join((
        "Number of MS1 scans: %d" %number_spectra_level1, 
        "Number of MS2 scans: %d" %number_spectra_level2,
        "Max retention time: %.2f sec" %max_rtime, 
    ))
    # props = dict(boxstyle='round', facecolor='y', alpha=0.1)
    ax[1].text(0.05, 0.95, note, transform=ax[1].transAxes, fontsize=8,
        verticalalignment='top',) # bbox=props)
    ax[1].set_xlabel("Retention time (sec)")
    ax[1].set_ylabel("Base peak intensity")
    ax[1].spines['top'].set_visible(False)
    ax[1].spines['right'].set_visible(False)
    ax[1].set_title("Ion chromatogram", fontsize=12)

    plt.tight_layout(pad=2.0)
    if outfile:
        plt.savefig(outfile)
        
        
        
if __name__ == "__main__":
    '''
    Example usage:
    python3 peek_ms_datafile.py -i B01_ss_017.mzML --title cool_data --output B01_ss_017.pdf
    '''
    
    import argparse
    parser = argparse.ArgumentParser(description='Plot summary of mass spec file')
    parser.add_argument('-i', type=str, help='Input file path')
    parser.add_argument('--title', type=str, default='', help='Figure title')
    parser.add_argument('--output', type=str, default="msfile_summary.pdf", help='Output file name.')
    args = parser.parse_args()

    summarize_ms_file(
        infile=args.i, 
        title=args.title,
        outfile=args.output
    )
