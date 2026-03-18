'''
Bottom part from ChatGPT, yet to validate. Using top SL simple function for now.

Common MSP Variations

| Variation             | Example                           |
| --------------------- | --------------------------------- |
| annotations           | `100 45 "loss of H2O"`            |
| comma-separated peaks | `100 45, 120 32, 140 10`          |
| missing `Num Peaks`   | peaks until blank line            |
| extra metadata        | `InChIKey`, `CAS`, `Comment`, etc |

'''
import re

#
# SL simple function
#

def parse_msp_to_listdict(file, field_separator=': ', return_peaks=True):
    '''
    parse a msp file within reasonable size.
    return list of dictionaries, [{key, value pairs}, ...]
    entry['peaks'] contains list of tuples if data are imported correctly. 
    '''
    features = []
    w = open(file).read().rstrip().split('\n\n')
    for block in w:
        d = {}
        lines = block.splitlines()
        data = []
        for line in lines:
            if field_separator in line:
                x = line.split(field_separator)
                d[x[0]] = x[1]
            elif line.strip():
                data.append(line)
        if return_peaks:
            try:
                peaks = []
                for line in data:
                    for x in line.split(';'):
                        entry = x.strip().split()
                        if len(entry) == 2:
                            peaks.append((float(entry[0]), float(entry[1])))
                d['peaks'] = peaks  
            except(TypeError, ValueError):
                print("Failed to convert peaks: ", entry, "Data: ", data)
        features.append(d)
        
    return features

#
# Below is ChatGPT code, yet to verify
#

def parse_peak_token(token):
    """
    Parse a single peak token like:
    '100:45'
    '100 45'
    '100:45:annotation'
    '100:45:"loss of H2O"'
    """
    parts = token.split(":")

    if len(parts) < 2:
        return None

    try:
        mz = float(parts[0])
        intensity = float(parts[1])
    except ValueError:
        return None

    annotation = ":".join(parts[2:]) if len(parts) > 2 else None

    return {
        "mz": mz,
        "intensity": intensity,
        "annotation": annotation
    }


def parse_peak_line(line):
    """
    Parse a line containing one or more peaks.
    Handles formats like:

    100 45
    100 45 "annotation"
    100:45
    100:45 120:30
    100 45, 120 30, 140 10
    """

    peaks = []

    # Remove commas separating peaks
    line = line.replace(",", " ")

    tokens = line.split()

    i = 0
    while i < len(tokens):

        # Case: mz intensity pair
        try:
            mz = float(tokens[i])
            intensity = float(tokens[i + 1])

            annotation = None

            if i + 2 < len(tokens) and not re.match(r'^-?\d+(\.\d+)?$', tokens[i + 2]):
                annotation = tokens[i + 2]
                i += 3
            else:
                i += 2

            peaks.append({
                "mz": mz,
                "intensity": intensity,
                "annotation": annotation
            })

        except (ValueError, IndexError):

            # Try colon style peak
            peak = parse_peak_token(tokens[i])
            if peak:
                peaks.append(peak)

            i += 1

    return peaks


def parse_msp(filepath):
    spectra = []
    current = None
    expected_peaks = None
    peak_count = 0

    with open(filepath, "r", encoding="utf-8", errors="ignore") as f:

        for raw_line in f:
            line = raw_line.strip()

            if not line:
                continue

            # metadata line
            if ":" in line and not re.match(r'^\d', line):

                key, value = line.split(":", 1)
                key = key.strip()
                value = value.strip()

                if key.lower() == "name":
                    if current:
                        spectra.append(current)

                    current = {"Name": value, "peaks": []}
                    expected_peaks = None
                    peak_count = 0
                    continue

                if current is None:
                    current = {"peaks": []}

                if key.lower() == "num peaks":
                    try:
                        expected_peaks = int(value)
                    except ValueError:
                        expected_peaks = None
                else:
                    current[key] = value

                continue

            # peak line
            if current is None:
                continue

            peaks = parse_peak_line(line)

            current["peaks"].extend(peaks)
            peak_count += len(peaks)

            if expected_peaks and peak_count >= expected_peaks:
                spectra.append(current)
                current = None
                expected_peaks = None
                peak_count = 0

    if current:
        spectra.append(current)

    return spectra


