import numpy as np
import pandas as pd
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import pymzml
import matplotlib.pyplot as plt
from collections import defaultdict
import seaborn as sns
import plotly.graph_objs as go
import os
import json

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

def generate_qc_report(job):
    mzml_file, output_file, spikeins = job

    if spikeins is None:
        print(f"Using default spike-ins for {mzml_file}")
        spikeins = [
            ('13C6-D-glucos', 187.0908, 0),
            ('trimethyl-13C3-caffeine', 198.0977, 0),
            ('15N-13C5-methionine', 156.0721, 0),
            ('13C5-L-glutamate', 153.0722, 0),
            ('15N2-uracil', 115.0286, 0),
            ('15N-L-tyrosine', 183.0782, 0),
        ]
    elif spikeins.endswith("json"):
        with open(spikeins) as f:
            spikeins = json.load(f)
    else:
        raise Exception("Spike-ins must be a list of tuples or a JSON file.")

    """
    Generates a QC report for an mzML file and saves it as an HTML file in the same location.
    """
    description_text = """
    The single file QC is a quick screening tool to check issues in mzML files.
    This includes: the TIC, histograms of intensity and mz, and examination of spike-in standards.
    """
    
    scan_levels = defaultdict(int)
    ions = defaultdict(int)
    exp = pymzml.run.Reader(mzml_file)
    
    for spec in exp:
        scan_levels[spec.ms_level] += 1
        if spec['positive_mode']:
            ions["Positive"] += 1
        else:
            ions["Negative"] += 1
    
    description_text += "Scan summary: " + ",".join([f"{v} scans for ms_level={k}" for k, v in scan_levels.items()]) + "\n"
    description_text += "Mode summary: " + ",".join([f"{v} scans for ion mode={k}" for k, v in ions.items()]) + "\n"
        
    def calcTIC(exp, mslevel=1):
        rt, tic = [], []
        scan_no = 0
        for spec in exp:
            if spec.ms_level == mslevel:
                scan_no += 1
                total = sum(int(peak[1]) for peak in spec.peaks('centroided'))
                tic.append(total)
                rt.append(scan_no)
        return rt, tic
    
    def extract_trio(exp, mslevel=1, min_intensity=1000):
        return [(int(peak[1]), peak[0], float(spec.scan_time[0]) * 60)
                for spec in exp if spec.ms_level == mslevel
                for peak in spec.peaks('centroided') if int(peak[1]) > min_intensity]
    
    def find_targets(trio_list, spikeins, mz_error=0.003, min_intensity=10000):
        solutions = {}
        for name, mz, _ in spikeins:
            solutions[name] = [trio for trio in trio_list if trio[0] > min_intensity and mz - mz_error < trio[1] < mz + mz_error]
        return solutions
    
    RTs, TICs = calcTIC(exp)
    df = pd.DataFrame({'x': RTs, 'y': TICs})
    fig_TIC = px.scatter(df, x='x', y='y', title="Total Ion Chromatogram")
    for_HTML = [fig_TIC]
    
    trio_list = extract_trio(exp)
    solutions = find_targets(trio_list, spikeins)
    summary_table = []
    
    for s, trios in solutions.items():
        df = pd.DataFrame({'mz': [x[1] for x in trios], 'intensity': [x[0] for x in trios], 'scan_no': [x[2] for x in trios]})
        target = next(sp for sp in spikeins if sp[0] == s)
        df['delta_mz'] = df['mz'] - target[1]
        
        fig = make_subplots(rows=1, cols=2, subplot_titles=['Mass Accuracy (Da)', 'TIC'])
        fig.add_trace(go.Scatter(x=df['delta_mz'], y=df['intensity'], mode='markers'), row=1, col=1)
        fig.add_trace(go.Scatter(x=df['scan_no'], y=df['intensity'], mode='markers'), row=1, col=2)
        
        max_intensity = max(trios, key=lambda x: x[0], default=(0, None))
        ppm_error = sum(df['delta_mz']) / len(df['delta_mz']) / target[1] * 1e6 if len(df) > 0 else "NA"
        detected = max_intensity[1] is not None
        summary_table.append({"target": s, "max_intensity": max_intensity[0] if detected else "NA",
                            "max_intensity_time": max_intensity[1] if detected else "NA", "ppm_error": ppm_error, "detected": detected})
        
        title = f"{s}<br>"
        if detected:
            title += f"max intensity of {max_intensity[0]} at rtime={max_intensity[1]}<br>"
            title += f"mean mz error of {ppm_error} ppm<br>"
        else:
            title += "NOT DETECTED!!!"
        fig.update_layout(title={'text': title, 'x': 0.5, 'xanchor': 'center'})
        for_HTML.append(fig)
    
    summary_df = pd.DataFrame(summary_table)
    html_table = summary_df.to_html()
    
    plots_html = "<div id='plots'>" + "".join(x.to_html(full_html=False) for x in for_HTML) + "</div>"
    formatted_text = description_text.replace('\n', '<br>')
    
    html_content = f"""
    <!DOCTYPE html>
    <html lang="en">
    <head>
        <meta charset="UTF-8">
        <title>Single File QC for {mzml_file}</title>
        <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    </head>
    <body>
        <h1>Single File QC for {mzml_file}</h1>
        <p>{formatted_text}</p>
        <div id="table">{html_table}</div>
        {plots_html}
    </body>
    </html>
    """
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    with open(output_file, 'w') as file:
        file.write(html_content)
    return output_file
    
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
    df = pd.read_table(infile, header=header, index_col=index_col, sep=sep)
    if len(df.columns) > max_col:
        print(f"[asari.qc] Warning: too many columns ({len(df.columns)}), truncating to {max_col}.")
        df = df.iloc[:, :max_col]
    return df