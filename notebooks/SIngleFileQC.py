import os
import plotly.express as px
from plotly.subplots import make_subplots
import pandas as pd
import sys
import pymzml
import matplotlib.pyplot as plt
from collections import defaultdict

# HTML content

description_text = """
The single file QC is a quick screening tool to check issues in mzML files.
This includes: the TIC, histograms of intensity and mz, and examination of spikein standards.
"""


scan_levels = defaultdict(int)
ions = defaultdict(int)
for spec in pymzml.run.Reader(sys.argv[1]):
    scan_levels[spec.ms_level] += 1
    if spec['positive_mode']:
        ions["Positive"] += 1
    else:
        ions["Negative"]
description_text += "Scan summary: " + ",".join([str(v) + " scans for ms_level=" + str(k) for k,v in scan_levels.items()]) + "\n"
description_text += "Mode summary: " + ",".join([str(v) + " scans for ion mode=" + str(k) for k,v in scan_levels.items()]) + "\n"
        
#
# key functions
#
# internal controls

# will change to using a tsv file of spike-in standards; and add negative mode

spikeins = [
    # name, M+H, RT,
    ('13C6-D-glucos', 187.0908, 0),
    ('trimethyl-13C3-caffeine', 198.0977, 0),
    ('15N-13C5-methionine', 156.0721, 0),
    ('13C5-L-glutamate', 153.0722, 0),
    ('15N2-uracil', 115.0286, 0),
    ('15N-L-tyrosine', 183.0782, 0),
]

other_controls = [
    # name, M+H, RT,
    ('phenylalanine', 166.0864, 48),
]



# other regularly detected endogeneous metabolites to be used for control
other_controls = [
    # name, M+H, RT,
    ('phenylalanine', 166.0864, 48),
]

def calcTIC(exp, mslevel=1):
    rt, tic = [], []
    scan_no = 0
    # Iterate through all spectra of the experiment
    for spec in exp:
        # Only calculate TIC for matching (MS1) spectra
        if spec.ms_level == mslevel:
            scan_no += 1
            total = 0
            for peak in spec.peaks('centroided'):
                total += int(peak[1])
            tic.append(total)
            rt.append(scan_no)
        else:
            scan_no += 1
    return rt, tic

def extract_trio(exp, mslevel=1, min_intensity=1000):
    trio_list = []
    for spec in exp:
        if spec.ms_level == mslevel:
            this_rt = spec.scan_time
            for peak in spec.peaks('centroided'):
                if int(peak[1]) > min_intensity:
                    trio_list.append((int(peak[1]), peak[0], float(this_rt[0]) * 60))
    return trio_list

        


for_HTML = []
exp = pymzml.run.Reader(sys.argv[1])

# plot TIC
RTs, TICs = calcTIC(exp)
df = pd.DataFrame({
    'x': RTs,
    'y': TICs,

})
fig_TIC = px.scatter(df, x='x', y='y', title="Total Ion Chromatogram")
for_HTML.append(fig_TIC)

# make spike-in plots


def find_targets(trio_list, spikeins, mz_error=0.003, min_intensity=10000):
    solutions = {}
    for spikein in spikeins:
        name, mz, rtime = spikein
        range_low, range_high = mz - mz_error, mz + mz_error
        good_trios = []
        for trio in trio_list:
            if trio[0] > min_intensity:
                if range_low < trio[1] < range_high:
                    good_trios.append(trio)
        solutions[name] = good_trios
    return solutions

trio_list = extract_trio(exp)
solutions = find_targets(trio_list, spikeins)
import plotly.graph_objs as go
summary_table = []
for s in solutions:
    df = pd.DataFrame({
        'mz': [x[1] for x in solutions[s]],
        'intensity': [x[0] for x in solutions[s]],
        'scan_no': [x[2] for x in solutions[s]]
    })
    target = [spikein for spikein in spikeins if spikein[0] == s][0]
    df['delta_mz'] = df['mz'] - target[1]


    fig = make_subplots(rows=1, cols=2, subplot_titles=['Mass Accuracy (Da)', 'TIC'])
    fig.add_trace(go.Scatter(x=df['delta_mz'], y=df['intensity'], mode='markers'), row=1, col=1)
    fig.add_trace(go.Scatter(x=df['scan_no'], y=df['intensity'], mode='markers'), row=1, col=2)
    
    title = s + "<br>"
    max_intensity = (0, None)
    for trio in solutions[s]:
        if trio[0] > max_intensity[0]:
            max_intensity = (trio[0], trio[2])
    if max_intensity[1]:
        title += "max intensity of " + str(max_intensity[0]) + " at rtime=" + str(max_intensity[1]) + "<br>"
        title += "mean mz error of " + str(sum(((df['delta_mz'])) / len(df['delta_mz']) / target[1]) * 1e6) + " ppm<br>"
        summary_table.append({"target": s, "max_intensity": max_intensity[0], "max_intensity_time": max_intensity[1], "ppm_error": sum(((df['delta_mz'])) / len(df['delta_mz']) / target[1]) * 1e6, "detected": True})
    else:
        title += "NOT DETECTED!!!"
        summary_table.append({"target": s, "max_intensity": "NA", "max_intensity_time": "NA", "ppm_error": "NA", "detected": False})
    fig.update_layout(title={'text': title, 'x': 0.5, 'xanchor': 'center'})

    fig.update_xaxes(title_text="m/z error (da)", row=1, col=1)
    fig.update_yaxes(title_text="intensity (absolute)", row=1, col=1)

    # Update axes labels for the second scatter plot
    fig.update_xaxes(title_text="scan time", row=1, col=2)
    fig.update_yaxes(title_text="intensity (absolute)", row=1, col=2)
    
    for_HTML.append(fig)

def highlight_ppm_error(val):
    if val == "NA":
        color = 'red'
    else:
        color = 'red' if float(val) < -5 or float(val) > 5 else 'black'
    return f'color: {color}'

def highlight_detected(val):
    if val:
        color = 'black'
    else:
        color = 'red'
    return f'color: {color}'

summary_df = pd.DataFrame(summary_table)
styled_df = summary_df.style.applymap(highlight_ppm_error, subset=['ppm_error']).applymap(highlight_detected, subset=['detected'])
html_table = styled_df.to_html()



# Sample data for the scatter plot


# Insert plot HTML into the main HTML content
plots_html = f"""
    <div id="plots">
""" 
for x in for_HTML:
    plots_html += "\t\t" + x.to_html(full_html=False)

plots_html += """
    </div>
"""

formatted_text = description_text.replace('\n', '<br>')


html_content = f"""
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Single File QC for $SAMPLE</title>
    <script src="https://cdn.plot.ly/plotly-latest.min.js"></script>
    <style>
        h1 {{
            text-align: center;
        }}
        #table-container {{
            display: flex;
            justify-content: center;
        }}
        table {{
            border-collapse: collapse;
            width: 50%;
        }}
        th, td {{
            border: 1px solid black;
            padding: 8px;
            text-align: left;
        }}
    </style>
</head>
<body>
    <h1>Single File QC for {sys.argv[1]}</h1>
    <p>{formatted_text}.</p>
    <div>
        <div id="table-container">
        <div id="table">
            {html_table}
        </div>
    </div>
    <div id="images">
        <!-- Images will go here -->
    </div>
    <div id="plots">
        <!-- Interactive plots will go here -->
    </div>
        <div id="json_data", style="display:none;">
        {pd.DataFrame(summary_table).to_json(orient="records")}
    </div>
</body>
</html>
"""

html_content = html_content.replace('<!-- Interactive plots will go here -->', plots_html)

# Save the HTML content to a file
output = os.path.abspath(sys.argv[1]).replace(".mzML", "_report.html")
with open(output, 'w') as file:
    file.write(html_content)

print("HTML file created successfully.")