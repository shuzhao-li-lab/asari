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

hv.extension('bokeh')
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
    '''
    Convert peak dictionary into readable HTML.
    May need to improve error handling since KeyError is potential problem.
    '''
    s = title
    info = [ ('id_number: ', d['id_number'], ' - ', 'parent_masstrack_id: ', d['parent_masstrack_id'], ' - ', 'parent_epd_id: ', d.get('parent_epd_id', '')),
             ('snr: ', d['snr'], ' - ', 'peak shape: ', round(d['goodness_fitting'],2), ' - ', 'cSelectivity: ', round(d['cSelectivity'],2),),
             ('height: ', d['height'], ' - ', 'peak_area: ', d['peak_area']),
             ('mz: ', round(d['mz'],4), ' - ', 'rtime: ', (round(d['rtime_left_base'],2), round(d['rtime'],2), round(d['rtime_right_base'],2))),
    ]
    for x in info:
        s += "<ul>" + ' '.join([str(ii) for ii in x]) + "</ul>"
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
# Summary panel
#
def get_summary_panel(project_desc, peakDict, epdDict, Ftable):
    desc0 = "Project retrieved from %s, %d features and %d empirical compounds." %(project_desc['outdir'],
                    len(peakDict), len(epdDict)
                    )
    html = "<h2>Data summary</h2><ul>" + desc0 + \
            "</ul><ul>This Dashboard shows all features of varying measurement quality. We recommend to use the preferred feature table for data analysis.</ul>"
    
    description = pn.pane.HTML(html)
    height = 240
    num_bins = int(np.sqrt(Ftable.shape[0]))
    # Feature m/z distribution histogram
    mz_distribution = Ftable['mz'].hvplot.hist(bins=num_bins, tools=[], 
                        height=height, alpha=0.7, color='gold', hover=False,
                        ylabel="Count", title="Feature distribution by m/z"
                        ).opts(toolbar=None)
    # eature distribution by RT scatterplot
    Ftable['peak_area_sqrt'] = np.sqrt(Ftable['peak_area'])
    rt_distribution = Ftable.hvplot.scatter(x='rtime', y='peak_area_sqrt', tools=[], height=height, size=2, hover=False,
                        title="Feature distribution by RT").opts(toolbar=None)
    # quality metrics
    Ftable['log10(SNR)'] = np.log10(Ftable['snr']+1)
    SNR = Ftable['log10(SNR)'].hvplot.hist(bins=num_bins, tools=[], 
                        height=height, alpha=0.7, color='green', hover=False, 
                        ylabel="Count", title="Signal to noise ratio"
                        ).opts(toolbar=None)
    PeakShape = Ftable['goodness_fitting'].hvplot.hist(bins=num_bins, tools=[], 
                        height=height, alpha=0.7, color='grey', hover=False, 
                        ylabel="Count", title="Peak shape"
                        ).opts(toolbar=None)
    cSel = Ftable['cSelectivity'].hvplot.hist(bins=num_bins, tools=[], 
                        height=height, alpha=0.7, color='red', hover=False, 
                        ylabel="Count", title="Chromatographic Selectivity"
                        ).opts(toolbar=None)

    feature_distribution = pn.Tabs(
        ("Feature m/z distribution histogram", mz_distribution),
        ("Feature distribution by RT scatterplot", rt_distribution),
        ("Signal to noise ratio", SNR),
        ("Peak shape", PeakShape),
        ("Chromatographic Selectivity", cSel),
        )
    return pn.Column(
                    description,
                    feature_distribution,
                    )


#
# Dashboard
#
def dashboard(project_desc, cmap, epd, Ftable):
    '''
    Panel based Dashboard.
    '''
    print("//*Asari dashboard*//   Press Control-C to exit.")
    peakDict, epdDict = epd_convert(epd)
    a_good_peak = find_a_good_peak(peakDict)
    track2peaks = track_to_peaks(peakDict)
    rt_list = [cmap['dict_scan_rtime'][ii] for ii in range(cmap['rt_length'])]
    summary = get_summary_panel(project_desc, peakDict, epdDict, Ftable)

    #
    # Feature browser
    #
    default_feature_number = int(a_good_peak['id_number'][1:])
    max_feature_number = len(peakDict) # asari feature ID starts with F1
    feature_num_slider = pn.widgets.IntSlider(name='Feature number', value=default_feature_number, 
                        start=1, end=max_feature_number, step=1,)
    feature_selector = pn.widgets.IntInput(name='', value=default_feature_number, 
                        start=1, end=max_feature_number) 
    feature_selector.link(feature_num_slider, value='value')
    feature_num_slider.link(feature_selector, value='value')

    def feature_info_by_feature_id(feature_number):
        try:
            info = convert_dict_html( peakDict[ 'F' + str(feature_number) ] )
        except KeyError:
            info = "<p>Feature info not found - %d.</p>" %feature_number
        return pn.pane.HTML(info)

    def cmapplot_track_by_feature_id(feature_number):
        id_number = 'F' + str(feature_number)
        _p = { 'height': 1000000, 'rtime_left_base': 0, 'rtime_right_base': 10, }
        p = peakDict.get(id_number, _p)
        height = max(1.3 * p['height'], 1000000)
        track_id_number = peakDict[id_number]['parent_masstrack_id']
        peak_Area = hv.Area(
            ([p['rtime_left_base'], p['rtime_right_base']], [height, height])).opts(
                color='red', alpha=0.1,
            )
        return cmapplot_mass_tracks(cmap, rt_list, 'blue', track_id_number) * peak_Area

        
    display_feature_info = pn.bind(
        feature_info_by_feature_id, feature_number=feature_selector,
        )
    plot_masstracks_fid = pn.bind(
        cmapplot_track_by_feature_id, feature_number=feature_selector,
        )
    
    feature_browser = pn.Column(
        "## Feature browser",
        pn.Row( pn.Column(feature_num_slider, feature_selector),
                display_feature_info ),
        plot_masstracks_fid,
        )

    #
    # m/z browser
    #
    default_mz = round(a_good_peak['mz'], 4)
    min_mz, max_mz = Ftable['mz'].min(), Ftable['mz'].max()
    mz_slider = pn.widgets.FloatSlider(name='Find by closest m/z', value=default_mz, 
                        step=0.001, start=min_mz, end=max_mz,)
    mz_selector = pn.widgets.FloatInput(name='', value=default_mz, 
                        step=0.001, start=min_mz, end=max_mz,)
    mz_selector.link(mz_slider, value='value')
    mz_slider.link(mz_selector, value='value')

    def get_features_by_mz(mz):
        # return tid_number, list of dictionary of features associated with the mass track
        tid_number = str(find_track_by_mz(cmap, rt_list, mz))
        try:
            return tid_number, [peakDict[F] for F in track2peaks[tid_number]]  
        except KeyError:
            return tid_number, []

    def track_info_by_mz(mz):
        tid_number, list_features = get_features_by_mz(mz)
        if list_features:
            info = ""
            for peak in list_features:
                info += "<p>" + convert_dict_html(peak, peak['id_number']) + "</p>"
        else:
            info = "<p>No qualified feature found on this mass track - %s.</p>" %tid_number
        return pn.pane.HTML(info)

    def cmapplot_track_by_mz(mz):
        tid_number, list_features = get_features_by_mz(mz)
        trackplot = cmapplot_mass_tracks(cmap, rt_list, 'green', tid_number)
        for p in list_features:
            height = max(1.3 * p['height'], 1000000)
            center = p['rtime']
            peak_Area = hv.Area(([center-1, center+1], [height, height])).opts(color='red', alpha=0.1)
            trackplot *= peak_Area
        return trackplot

    display_track_info = pn.bind(
        track_info_by_mz, mz=mz_selector,
        )
    plot_masstracks_mz = pn.bind(
        cmapplot_track_by_mz, mz=mz_selector,
        )

    mz_browser = pn.Column(
        "## View composite mass track by m/z",
        pn.Row( pn.Column(mz_slider, mz_selector),
                display_track_info ),
        plot_masstracks_mz,
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
