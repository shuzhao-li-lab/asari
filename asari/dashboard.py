"""
Functions for subcommand `viz`
"""
import os
import html
import bisect
import pickle
import json
import numpy as np
import pandas as pd

import hvplot.pandas
import panel as pn
import holoviews as hv

hv.extension('bokeh')
pn.extension(sizing_mode="stretch_width")


FEATURE_INFO_COLUMNS = ['id_number', 'mz', 'rtime', 'rtime_left_base', 'rtime_right_base',
                        'parent_masstrack_id', 'peak_area', 'cSelectivity', 'goodness_fitting',
                        'snr', 'detection_counts', 'RI']


def get_data_start_col(columns):
    """Index of first sample column in an asari feature table: after 'RI' if present, else after 'detection_counts'."""
    columns = list(columns)
    if 'RI' in columns:
        return columns.index('RI') + 1
    return columns.index('detection_counts') + 1


def read_feature_table(infile, load_sample_limit=20):
    """
    Read an asari feature table (tsv) using id_number as index.
    Sample columns are truncated to the first load_sample_limit.
    """
    header = pd.read_csv(infile, sep='\t', nrows=0).columns
    usecols = list(header[: get_data_start_col(header) + load_sample_limit])
    return pd.read_csv(infile, sep='\t', usecols=usecols, index_col=0, header=0)


def read_project(datadir, load_sample_limit=20):
    """
    Get all project data. Requires a project processed with `--diagnosis True`, which exports cmap.pickle.

    Returns
    -------
    project_desc :
        dict of project meta data
    cmap :
        composite map, ['_number_of_samples_', 'rt_length', 'rt_reference_landmarks', 'rt_records',
        'dict_scan_rtime', 'list_mass_tracks', 'MassGrid']
    Ftable :
        pandas dataframe of full feature table, id_number as index. Truncated if samples more than load_sample_limit.
    Ptable :
        pandas dataframe of preferred feature table, same format as Ftable.
    """
    datadir = os.path.abspath(datadir)
    project_desc = json.load(open(os.path.join(datadir, 'project.json')))
    cmap_file = os.path.join(datadir, 'export', 'cmap.pickle')
    if not os.path.exists(cmap_file):
        raise FileNotFoundError(f"{cmap_file} not found. Run `asari process` with `--diagnosis True` to export data for viz.")
    cmap = pickle.load(open(cmap_file, 'rb'))
    Ftable = read_feature_table(os.path.join(datadir, 'export', 'full_Feature_table.tsv'), load_sample_limit)
    Ptable = read_feature_table(os.path.join(datadir, 'preferred_Feature_table.tsv'), load_sample_limit)
    return project_desc, cmap, Ftable, Ptable


def table_to_feature_dict(Ftable, cmap, rt_list):
    """
    Convert feature table to {id_number: feature dict}, using feature columns only.
    'height' is added as max composite intensity of the parent mass track within the peak RT range.
    """
    cols = [c for c in FEATURE_INFO_COLUMNS[1:] if c in Ftable.columns]
    rt = np.asarray(rt_list)
    peakDict = {}
    for fid, F in Ftable[cols].to_dict('index').items():
        F['id_number'] = fid
        track = cmap['list_mass_tracks'].get(int(F['parent_masstrack_id']))
        height = 0
        if track is not None:
            left = np.searchsorted(rt, F['rtime_left_base'], side='left')
            right = np.searchsorted(rt, F['rtime_right_base'], side='right')
            segment = np.asarray(track['intensity'])[left:max(right, left + 1)]
            height = int(segment.max()) if len(segment) else 0
        F['height'] = height
        peakDict[fid] = F
    return peakDict


def plot_xic(xics, mz_dict, track_id):
    """
    Generates scatter plot for a mass track as dataframe. Test function.
    """
    title = "Mass track viewer - m/z %4.4f" %(mz_dict[track_id])
    plot = xics.hvplot.scatter(x='rt', y=track_id, height=400, title=title,)
    plot.opts(active_tools=['box_zoom'], ylabel="Normalized intensity", )
    return plot

def cmapplot_mass_tracks(cmap, rt_list, color, track_id_number):
    """
    return a hv plot of mass track, by track_id_number

    Parameters
    ----------
    cmap :
        CMAP as from ext_Experiment {
            '_number_of_samples_': self.CMAP._number_of_samples_,
            'rt_length': self.CMAP.rt_length,
            'dict_scan_rtime': self.CMAP.dict_scan_rtime,
            'list_mass_tracks': self.CMAP.composite_mass_tracks,
            'MassGrid': dict(self.CMAP.MassGrid),
        }
    rt_list :
        list of retention time
    color :
        color for scatter
    track_id_number :
        index for a mass track in cmap['list_mass_tracks']

    Returns
    -------
    A holoviz plot object.
    """
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
    """
    Convert a peak dictionary into readable HTML.
    May need to improve error handling since KeyError is potential problem.
    Returns HTML as a string.
    """
    s = title
    info = [ ('id_number: ', d['id_number'], ' - ', 'parent_masstrack_id: ', d['parent_masstrack_id'], ' - ',
              'detection_counts: ', d.get('detection_counts', '')),
             ('snr: ', d['snr'], ' - ', 'peak shape: ', round(d['goodness_fitting'],2), ' - ', 'cSelectivity: ', round(d['cSelectivity'],2),),
             ('height: ', d['height'], ' - ', 'peak_area: ', d['peak_area']),
             ('mz: ', round(d['mz'],4), ' - ', 'rtime: ', (round(d['rtime_left_base'],2), round(d['rtime'],2), round(d['rtime_right_base'],2))),
    ]
    if d.get('RI') is not None and not pd.isna(d.get('RI')):
        info.append(('RI: ', d['RI'],))
    for x in info:
        s += "<ul>" + ' '.join([str(ii) for ii in x]) + "</ul>"
    return s


def convert_dict_markdown(d, title=''):
    """
    Convert a peak dictionary into Markdown string.
    """
    s = title
    for k,v in d.items():
        if k not in ['apex', 'left_base', 'right_base']:
            s += '- ' + k + ': \t' + str(v) + '\n'
    return s + '\n'

def track_to_peaks(peakDict):
    """
    Parameters
    -----
    peakDict, peak dictionary using 'id_number' as keys.

    Returns
    -------
    a dictionary using parent_masstrack_id as keys and ['id_number'] as values.
    """
    t = {}
    for P in peakDict.values():
        mid = int(P["parent_masstrack_id"])
        if mid in t:
            t[mid].append(P["id_number"])
        else:
            t[mid] = [P["id_number"]]
    return t

def find_track_by_mz(cmap, rt_list, mz):
    """
    return track_id_number by cloesest m/z.
    """
    L = [(abs(mz-T['mz']), T['id_number']) for T in cmap['list_mass_tracks'].values()]
    L2 = [x for x in L if x[0] < 0.1]
    if L2:
        return sorted(L2)[0][1]
    else:
        return sorted(L)[0][1]

def find_a_good_peak(peakDict):
    """
    To find a good example feature/peak, which can be used as the feature on landing page.
    """
    good = [P for P in peakDict.values() if P['goodness_fitting'] > 0.9 and P['cSelectivity'] > 0.9]
    if good:
        return good[0]
    return max(peakDict.values(), key=lambda P: P['peak_area'])

def prepare_rt_alignment(cmap):
    """this method prepares necessary dataframe for retention time alignment.

    Parameters
    ----------
    cmap : dict
        pickle information read include all we need for retention time alignment

    Returns
    -------
    dataframe
        pandas dataframe ready for holoview drawing

    """
    sample_reference_list = [rt_record['reverse_rt_cal_dict'] for rt_record in cmap['rt_records']]
    rt_length = cmap['rt_length']

    for i, d in enumerate(sample_reference_list):
        # imputate the equal key-value pair
        missing_keys = {k: k for k in range(rt_length) if k not in d}
        d.update(missing_keys)

        # sort the dicts by key & map scan number to rtime & calculate the deviation
        deviation_dict = {cmap['dict_scan_rtime'][k]: cmap['dict_scan_rtime'][d[k]]-cmap['dict_scan_rtime'][k] for k in sorted(d)}

        sample_reference_list[i] = deviation_dict

    sample_name_list = [record['name'] for record in cmap['rt_records']]
    deviation_df = pd.DataFrame(sample_reference_list).T
    deviation_df.reset_index().rename(columns={'index':'rtime'})
    # rename columns
    deviation_df.columns = sample_name_list

    return deviation_df

#
# Summary panel
#
def get_summary_panel(project_desc, Ftable, cmap, table_name='full'):
    """
    Get a summary panel, returns a panel.Column of multiple tabs for summary metrics.
    """
    desc0 = "Project retrieved from %s, %d features in the %s feature table." %(project_desc['outdir'],
                    Ftable.shape[0], table_name
                    )
    html = "<h2>Data summary</h2><ul>" + desc0 + \
            "</ul><ul>This Dashboard shows all features of varying measurement quality. We recommend to use the preferred feature table for data analysis.</ul>"

    description = pn.pane.HTML(html)
    height = 330
    width = 1000

    #metadata = []
    #for key, value in project_desc.items():
    #    metadata.append({"Parameter": key, "Value": value})
    #project_metadata = pd.DataFrame(metadata)
    #metadata_table = project_metadata.hvplot.table(columns=['Parameter', 'Value'], width=width, height=height).opts(toolbar=None)

    num_bins = int(np.sqrt(Ftable.shape[0]))
    # each plot keeps its own axes (not linked across tabs) and the default toolbar (pan, zoom, reset, save)
    hist_opts = dict(bins=num_bins, width=width, height=height, alpha=0.7, hover=False, ylabel="Count")
    # Feature m/z distribution histogram
    mz_distribution = Ftable['mz'].hvplot.hist(color='gold', title="Feature distribution by m/z", **hist_opts)
    # Feature distribution by RT scatterplot
    Ftable['peak_area_sqrt'] = np.sqrt(Ftable['peak_area'])
    rt_distribution = Ftable.hvplot.scatter(x='rtime', y='peak_area_sqrt', width=width, height=height, size=2,
                        hover=False, title="Feature distribution by RT")
    # quality metrics
    Ftable['log10(SNR)'] = np.log10(Ftable['snr']+1)
    SNR = Ftable['log10(SNR)'].hvplot.hist(color='green', title="Signal to noise ratio", **hist_opts)
    PeakShape = Ftable['goodness_fitting'].hvplot.hist(color='grey', title="Peak shape", **hist_opts)
    cSel = Ftable['cSelectivity'].hvplot.hist(color='red', title="Chromatographic Selectivity", **hist_opts)

    # hard cap for legend is 20
    if_legend = cmap['_number_of_samples_'] <= 20
    RTAlign = prepare_rt_alignment(cmap).hvplot.line(title="Retention Time Deviation vs Retention Time",xlabel="Retention Time (sec)",
                ylabel="Retention Time Deviation (sec)", width=width, height=height,hover=True).opts(
                show_legend=if_legend, legend_opts=dict(
                title='',label_text_font_size='8pt',spacing=-5, location=(5, (height/2-cmap['_number_of_samples_']*11)), padding=3))

    def _pane(plot):
        return pn.pane.HoloViews(plot.opts(shared_axes=False), linked_axes=False)

    feature_distribution = pn.Tabs(
        ("Feature m/z distribution histogram", pn.Column(_pane(mz_distribution))),
        ("Feature distribution by RT scatterplot", pn.Column(_pane(rt_distribution))),
        ("Signal to noise ratio", pn.Column(_pane(SNR))),
        ("Peak shape", pn.Column(_pane(PeakShape))),
        ("Chromatographic Selectivity", pn.Column(_pane(cSel))),
        ("Retention Time Alignment", pn.Column(_pane(RTAlign))),
        )
    return pn.Column(
                    description,
                    feature_distribution,
                    )


#
# Dashboard
#
def build_dashboard(project_desc, cmap, Ftable, table_name='full'):
    """
    Build the panel based Dashboard. Returns a panel template.

    Features are retrieved from Ftable (full or preferred feature table, id_number as index);
    mass tracks from cmap.
    """
    rt_list = [cmap['dict_scan_rtime'][ii] for ii in range(cmap['rt_length'])]
    peakDict = table_to_feature_dict(Ftable, cmap, rt_list)
    a_good_peak = find_a_good_peak(peakDict)
    track2peaks = track_to_peaks(peakDict)
    summary = get_summary_panel(project_desc, Ftable, cmap, table_name)

    #
    # Feature browser
    #
    # feature numbers present in the table (asari feature ID is 'F' + number); not continuous in preferred table
    feature_numbers = sorted(int(fid[1:]) for fid in peakDict)
    default_feature_number = int(a_good_peak['id_number'][1:])
    feature_num_slider = pn.widgets.DiscreteSlider(name='Feature number', options=feature_numbers,
                        value=default_feature_number)
    feature_selector = pn.widgets.IntInput(name='', value=default_feature_number,
                        start=feature_numbers[0], end=feature_numbers[-1])

    def _available_feature_number(n, previous=None):
        """
        Next available feature number >= n (or the last one).
        If stepping down by 1 from previous (e.g. down arrow), the previous available number <= n instead.
        """
        if previous is not None and n == previous - 1:
            ii = bisect.bisect_right(feature_numbers, n)
            return feature_numbers[max(ii - 1, 0)]
        ii = bisect.bisect_left(feature_numbers, n)
        return feature_numbers[min(ii, len(feature_numbers) - 1)]

    def _on_selector(event):
        # move a number not in the table to the next available feature (previous one on a step down), then sync slider
        n = _available_feature_number(event.new, event.old)
        if n != event.new:
            feature_selector.value = n
        elif feature_num_slider.value != n:
            feature_num_slider.value = n

    def _on_slider(event):
        if feature_selector.value != event.new:
            feature_selector.value = event.new

    feature_selector.param.watch(_on_selector, 'value')
    feature_num_slider.param.watch(_on_slider, 'value')

    def feature_info_by_feature_id(feature_number):
        try:
            info = convert_dict_html( peakDict[ 'F' + str(feature_number) ] )
        except KeyError:
            info = "<p>Feature info not found - %d.</p>" %feature_number
        return pn.pane.HTML(info)

    def cmapplot_track_by_feature_id(feature_number):
        id_number = 'F' + str(feature_number)
        if id_number not in peakDict:
            return hv.Curve([]).opts(title="Feature not found - %s." %id_number, height=400)
        p = peakDict[id_number]
        height = max(1.3 * p['height'], 1000000)
        track_id_number = p['parent_masstrack_id']
        peak_Area = hv.Area(
            ([p['rtime_left_base'], p['rtime_right_base']], [height, height])).opts(
                color='red', alpha=0.1,
            )
        return cmapplot_mass_tracks(cmap, rt_list, 'blue', track_id_number) * peak_Area


    display_feature_info = pn.bind(
        feature_info_by_feature_id, feature_number=feature_num_slider,
        )
    plot_masstracks_fid = pn.bind(
        cmapplot_track_by_feature_id, feature_number=feature_num_slider,
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
        tid_number = int(find_track_by_mz(cmap, rt_list, mz))
        return tid_number, [peakDict[F] for F in track2peaks.get(tid_number, [])]

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

    #
    # Co-elution features
    #

    min_rt, max_rt = Ftable['rtime'].min(), Ftable['rtime'].max()
    min_mz, max_mz = Ftable['mz'].min(), Ftable['mz'].max()
    rt_range_slider = pn.widgets.FloatSlider(name='coelute_tol', value=1, step=0.05, start=0, end=10)
    rt_slider = pn.widgets.FloatSlider(name='rtime', value=min_rt * 1.05, step=0.01, start=min_rt, end=max_rt)

    def features_by_rt_with_lines(rtime, rtime_tol):
        # Filter features based on retention time
        Fsub = Ftable[Ftable['rtime'].between(rtime - rtime_tol, rtime + rtime_tol)]
        try:
            # Prepare the data for horizontal lines
            segments_data = pd.DataFrame({
                'x': Fsub['mz'],
                'y0': 0,                # Start Y position is always 0
                'y1': Fsub['peak_area'], # End Y position is the peak area
            })

            # Create a horizontal line plot using hv.Segments
            return hv.Segments(segments_data, kdims=['x', 'y0', 'x', 'y1']).opts(
                xlabel="m/z",
                ylabel="Peak Area",
                title=f"Features near RT: {rtime} ± {rtime_tol}",
                color="blue",
                line_width=2,
                active_tools=['box_zoom'],
                xlim=(min_mz * .95, max_mz * 1.05)
            )
        except:
            return hv.Curve([]).opts(title="No features found in this range.", height=400)

    display_frag = pn.bind(features_by_rt_with_lines, rtime=rt_slider, rtime_tol=rt_range_slider)    # Create the RT browser interface
    frag_browser = pn.Column(
        "## Co-elution features",
        pn.pane.HTML('''<p> Useful for interpretation of isotopologues, adducts and fragments etc. </p>
                    '''),
        pn.Row(rt_slider, rt_range_slider),
        display_frag,
    )

    disclaimer = pn.pane.HTML('''<p>Dangerous data. Use at your own risk.</p>
    <p>Asari source code is hosted at <a href="https://github.com/shuzhao-li/asari" target="_blank">https://github.com/shuzhao-li/asari</a>.
    Feel free to contribute, use GitHub Issues to report bugs and request features.</p>
    ''')

    params_text = '; '.join(f"{k}: {v}" for k, v in project_desc.items())
    project_params = pn.pane.HTML("<h2>Project Parameters</h2><p>" + html.escape(params_text) + "</p>")

    template = pn.template.FastListTemplate(
        site = "Asari Dashboard",
        title = project_desc["project_name"],
        main = [
                summary,
                feature_browser,
                mz_browser,
                frag_browser,
                project_params,
                disclaimer,
                ],
        )

    return template


def dashboard(project_desc, cmap, Ftable, table_name='full'):
    """
    Launch the panel based Dashboard in web browser.
    """
    print("//*Asari dashboard*//   Press Control-C to exit.")
    pn.serve(build_dashboard(project_desc, cmap, Ftable, table_name))


#
# -----------------------------------------------------------------------------
#

if __name__ == '__main__':
    # datadir = 'asari_project_MT01_51162751'

    import sys
    datadir = sys.argv[1]

    project_desc, cmap, Ftable, Ptable = read_project(datadir)
    dashboard(project_desc, cmap, Ftable)
