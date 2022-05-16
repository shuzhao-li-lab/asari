'''
Functions for subcommand `viz`
'''
import os
import pickle
import json
import numpy as np
import pandas as pd

import hvplot.pandas
import panel as pn
import holoviews as hv

pn.extension(sizing_mode="stretch_width")


def epd_convert(epd_dict):
    '''
    Format epd_dict to two dicts: peakDict, epdDict
    '''
    peakDict, epdDict = {}, {}
    for k,v in epd_dict.items():
        for P in v['MS1_pseudo_Spectra']:
            peakDict[P['id_number']] = P
        v['MS1_pseudo_Spectra'] = [P['id_number'] for P in v['MS1_pseudo_Spectra']]
        epdDict[k] = v
    return peakDict, epdDict

def read_project(datadir, load_sample_limit=20):
    '''
    Get all project data.
    Return
    ======
    project_desc: dict of project meta data
    cmap: composite map, ['_number_of_samples_', 'rt_length', 'dict_scan_rtime', 'list_mass_tracks', 'MassGrid']
    epd: dict of empirical compounds
    Ftable: pandas dataframe of feature table. Truncated if samples more than load_sample_limit.
    '''
    project_desc = json.load(open(os.path.join(datadir, 'project.json')))
    cmap = pickle.load( open(os.path.join(datadir, 'export', 'cmap.pickle'), 'rb') )
    # xics, mz_dict, massgrid = reformat_cmap()

    epd = pickle.load( open(os.path.join(datadir, 'export', 'epd.pickle'), 'rb') )
    if 'number_of_samples' in project_desc and project_desc['number_of_samples'] > load_sample_limit:
        Ftable = pd.read_csv( os.path.join(datadir, 'export', 'full_Feature_table.tsv'), 
                                sep='\t', index_col=0, header=0, usecols=range(10 + load_sample_limit) )
    else:
        Ftable = pd.read_csv( os.path.join(datadir, 'export', 'full_Feature_table.tsv'), 
                                sep='\t', index_col=0, header=0 )

    return project_desc, cmap, epd, Ftable

def plot_xic(xics, mz_dict, track_id):
    '''
    track_id is a str
    '''
    title = "Mass track viewer - m/z %4.4f" %(mz_dict[track_id])
    plot = xics.hvplot.scatter(x='rt', y=track_id, height=400, title=title,)
    plot.opts(active_tools=['box_zoom'], ylabel="Normalized intensity", )
    return plot

def cmapplot_mass_tracks(cmap, rt_list, color, track_id_number):
    '''
    return a hv plot of mass track, by track_id_number
    '''
    track_id_number = int(track_id_number)
    mz = cmap['list_mass_tracks'][track_id_number]['mz']
    title = "Mass track viewer - (track_id %d, m/z %4.4f)" %(track_id_number, mz)
    data = pd.DataFrame({
                    'rt': rt_list, 'intensity': cmap['list_mass_tracks'][track_id_number]['intensity']
                    })
    plot = data.hvplot.scatter(x='rt', y='intensity', color=color, height=400, title=title, size=8)
    plot.opts(active_tools=['box_zoom'], ylabel="Intensity (composite)", )
    return plot

def convert_dict_html(d, title=''):
    s = title
    for k,v in d.items():
        if k not in ['apex', 'left_base', 'right_base']:
            s += '<ul>' + k + ': \t' + str(v) + '</ul>'
    return s

def convert_dict_markdown(d, title=''):
    s = title
    for k,v in d.items():
        if k not in ['apex', 'left_base', 'right_base']:
            s += '- ' + k + ': \t' + str(v) + '\n'
    return s + '\n'

def track_to_peaks(peakDict):
    t = {}
    for P in peakDict.values():
        mid = str(P["parent_masstrack_id"])
        if mid in t:
            t[mid].append(P["id_number"])
        else:
            t[mid] = [P["id_number"]]
    return t

def find_track_by_mz(cmap, rt_list, mz):
    '''
    return track_id_number by cloesest m/z.
    '''
    L = [(abs(mz-T['mz']), T['id_number']) for T in cmap['list_mass_tracks'].values()]
    L2 = [x for x in L if x[0] < 0.1]
    if L2:
        return sorted(L2)[0][1]
    else:
        return sorted(L)[0][1]

def find_a_good_peak(peakDict):
    '''find a good example feature/peak'''
    good = [P for P in peakDict.values() if P['goodness_fitting'] > 0.9 and P['cSelectivity'] > 0.9]
    return good[0]

#
# Dashboard
#
def dashboard(project_desc, cmap, epd, Ftable):
    '''
    Panel based Dashboard
    '''
    print("//*Asari dashboard*//   Press Control-C to exit.")
    peakDict, epdDict = epd_convert(epd)
    a_good_peak = find_a_good_peak(peakDict)
    track2peaks = track_to_peaks(peakDict)
    rt_list = [cmap['dict_scan_rtime'][ii] for ii in range(cmap['rt_length'])]

    desc0 = "Project retrieved from %s, %d features and %d empirical compounds." %(project_desc['outdir'],
                    len(peakDict), len(epdDict)
                    )
    description = pn.pane.HTML("<h2>Data summary</h2><p>" + desc0 + "</p>")
    disclaimer = pn.pane.HTML('''<p>Dangerous data. Use at your own risk.</p>
    <p>Asari source code is hosted at <a href="https://github.com/shuzhao-li/asari" target="_blank">https://github.com/shuzhao-li/asari</a>.
    Feel free to contribute, use GitHub Issues to report bugs and request features.</p>
    ''')

    num_bins = int(np.sqrt(Ftable.shape[0]))
    mz_distribution = Ftable['mz'].hvplot.hist(bins=num_bins, tools=[], width=500, height=200, hover=False,
                        ylabel="Count",
                        title="Feature distribution by m/z")
    Ftable['peak_area_sqrt'] = np.sqrt(Ftable['peak_area'])
    rt_distribution = Ftable.hvplot.scatter(x='rtime', y='peak_area_sqrt', tools=[], width=500, height=200, size=2, hover=False,
                        title="Feature distribution by RT")
    feature_distribution = mz_distribution + rt_distribution
    feature_distribution.opts(toolbar=None)
    summary = pn.Column(
        description,
        feature_distribution,
        )

    #
    # Feature browser
    #
    default_feature_number = int(a_good_peak['id_number'][1:])
    feature_selector = pn.widgets.IntInput(name='Feature number', value=default_feature_number, 
                        start=1, end=len(peakDict)) # asari feature ID starts with F1

    def feature_info_by_feature_id(feature_number):
        return pn.pane.Markdown(convert_dict_markdown(peakDict[ 'F' + str(feature_number) ]))

    def cmapplot_track_by_feature_id(feature_number):
        color = 'blue'
        track_id_number = peakDict['F' + str(feature_number)]['parent_masstrack_id']
        return cmapplot_mass_tracks(cmap, rt_list, color, track_id_number)
        
    display_feature_info = pn.bind(
        feature_info_by_feature_id, feature_number=feature_selector,
        )
    plot_masstracks_fid = pn.bind(
        cmapplot_track_by_feature_id, feature_number=feature_selector,
        )
    
    feature_browser = pn.Column(
        "## Feature browser",
        pn.Row( feature_selector, display_feature_info ),
        plot_masstracks_fid,
        )

    #
    # m/z browser
    #
    default_mz = round(a_good_peak['mz'], 4)
    mz_selector = pn.widgets.FloatInput(name='Search m/z', value=default_mz, 
                        step=0.001, start=Ftable['mz'].min(), end=Ftable['mz'].max()
                        )
    
    def track_info_by_mz(mz):
        tid_number = str(find_track_by_mz(cmap, rt_list, mz))
        info = ""
        for F in track2peaks[tid_number]:
            info += "<p>" + F + " = " + str(peakDict[F]) + "</p>"
        return pn.pane.HTML(info)

    def cmapplot_track_by_mz(mz):
        color = 'green'
        track_id_number = str(find_track_by_mz(cmap, rt_list, mz))
        return cmapplot_mass_tracks(cmap, rt_list, color, track_id_number)
        
    display_track_info = pn.bind(
        track_info_by_mz, mz=mz_selector,
        )
    plot_masstracks_mz = pn.bind(
        cmapplot_track_by_mz, mz=mz_selector,
        )

    mz_browser = pn.Column(
        "## Mass track by m/z",
        pn.Row( mz_selector, display_track_info ),
        plot_masstracks_mz,
        )
    
    template = pn.template.FastListTemplate(
        site = "Asari Dashboard", 
        title = project_desc["project_name"],
        main = [
                summary,
                feature_browser,
                mz_browser,
                disclaimer,
                ], 
        )

    pn.serve(template)



#
# -----------------------------------------------------------------------------
#




if __name__ == '__main__':
    # datadir = 'asari_project_MT01_51162751'

    import sys
    datadir = sys.argv[1]
    
    project_desc, cmap, epd, Ftable = read_project(datadir)
    dashboard(project_desc, cmap, epd, Ftable)
