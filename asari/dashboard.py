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


# datadir = 'asari_project_MT01_51162751'

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


def cmapplot_mass_tracks(cmap, rt_list, track_id_number):
    track_id_number = int(track_id_number)
    mz = cmap['list_mass_tracks'][track_id_number]['mz']
    title = "Mass track viewer - (track_id %d, m/z %4.4f)" %(track_id_number, mz)
    data = pd.DataFrame({
                    'rt': rt_list, 'intensity': cmap['list_mass_tracks'][track_id_number]['intensity']
                    })
    plot = data.hvplot.scatter(x='rt', y='intensity', height=400, title=title,)
    plot.opts(active_tools=['box_zoom'], ylabel="Intensity (composite)", )
    return plot


def convert_dict_html(d, title):
    s = "<h3>%s</h3>" %title
    # for k in ['mz', 'rtime', 'rtime_left_base', 'rtime_right_base', ]
    for k,v in d.items():
        s += '<ul>' + k + ': \t' + str(v) + '</ul>'
    return s



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



def dashboard(project_desc, cmap, epd, Ftable):

    peakDict, epdDict = epd_convert(epd)
    a_good_peak = find_a_good_peak(peakDict)
    track2peaks = track_to_peaks(peakDict)

    desc0 = "Project retrieved from %s, %d features and %d empirical compounds." %(project_desc['outdir'],
                    len(peakDict), len(epdDict)
                    )
    description = pn.pane.HTML("<h3>Data summary</h3><p>" + desc0 + "</p>")

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

    # asari feature ID starts with F1
    feature_selector = pn.widgets.IntInput(name='Feature number', value=99, 
                        start=1, end=len(peakDict))
    mz_selector = pn.widgets.FloatInput(name='Search m/z', value=473.0565, step=0.001, 
                        start=Ftable['mz'].min(), end=Ftable['mz'].max()
                        )


    track_id_number = a_good_peak['parent_masstrack_id']

    rt_list = [cmap['dict_scan_rtime'][ii] for ii in range(cmap['rt_length'])]
    plot_masstracks = cmapplot_mass_tracks(cmap, rt_list, track_id_number=track_id_number)
    
    pstr = '\n'.join(
        [convert_dict_html(peakDict[x], x) for x in track2peaks[ str(track_id_number) ]
        ])
    feature_info = pn.pane.Markdown(pstr)


    data_block = pn.Column(
        pn.Row( mz_selector, feature_selector ),
        plot_masstracks,
        feature_info,
        )

    disclaimer = pn.pane.HTML('''<p>Dangerous data. Use at your own risk.</p>
    <p>Asari source code is hosted at <a href="https://github.com/shuzhao-li/asari" target="_blank">https://github.com/shuzhao-li/asari</a>.
    Feel free to contribute, use GitHub Issues to report bugs and request features.</p>
    ''')
    
    template = pn.template.FastListTemplate(
        site = "Asari Dashboard", 
        title = project_desc["project_name"],
        main = [
                summary,
                data_block,
                disclaimer
                ], 
        )


    pn.serve(template)



if __name__ == '__main__':
    print("//*Asari dashboard*//   Press Control-C to exit.")
    project_desc, cmap, epd, Ftable = read_project(datadir)
    dashboard(project_desc, cmap, epd, Ftable)

