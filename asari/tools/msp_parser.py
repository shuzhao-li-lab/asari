"""
MSP and MGF parsers.
"""
import re


# Controlled vocabulary: canonical key -> variants seen in .msp files.
# Keys are matched after normalization (case-insensitive, ignoring spaces, '_', '.', '-'),
# so 'Retention_index', 'RETENTION INDEX' and 'retention-index' are the same key.
# Order matters - the first variant found in an entry wins,
# so preferred variants go first (e.g. exact mass before MW, which is often nominal).
# 'peaks' is not listed as parse_msp_to_listdict already writes entry['peaks'].
# RTINSECONDS is deliberately left out to avoid unit mix-up.
# 'RI' (retention index) is kept separate from 'RETENTIONTIME' and converted to float;
# note that some GC libraries (e.g. gcms_corelib) store RI under RETENTIONTIME.
_MSP_variants = {
    'Name': ['Name', 'Compound Name', 'Compound_Name', 'CompoundName', 'Compound', 'Title'],
    'InChIKey': ['InChIKey', 'InchiKey', 'InChI Key', 'InChI_Key', 'Inchi_Key'],
    'ExactMass': ['ExactMass', 'Exact Mass', 'Exact_Mass', 'Monoisotopic Mass', 'Monoisotopic_Mass',
                  'MonoisotopicMass', 'Exact Mol Weight', 'MW', 'Mol Weight', 'Molecular Weight',
                  'Molecular_Weight', 'MolecularWeight'],
    'Formula': ['Formula', 'Molecular Formula', 'Molecular_Formula', 'MolecularFormula', 'Mol Formula'],
    'RETENTIONTIME': ['RETENTIONTIME', 'RetentionTime', 'Retention Time', 'Retention_Time',
                      'RT', 'RT_min', 'RTINMINUTES', 'Retention Time (min)'],
    'RI': ['RI', 'Retention_index', 'Retention Index', 'RetentionIndex', 'RETENTIONINDEX',
           'Kovats', 'Kovats_RI', 'Kovats RI', 'Kovats Index', 'Kovats_Index', 'KovatsIndex',
           'Kovats Retention Index', 'RI_Kovats', 'KI', 'LRI', 'Linear Retention Index'],
    'SMILES': ['SMILES', 'Canonical SMILES', 'Canonical_SMILES', 'Isomeric SMILES', 'Computed SMILES'],
    'PrecursorMZ': ['PrecursorMZ', 'Precursor_MZ', 'Precursor MZ', 'PEPMASS', 'Precursor m/z'],
    'Num Peaks': ['Num Peaks', 'Num_Peaks', 'NumPeaks', 'Num.peaks', 'Number of Peaks', 'NPeaks'],
}

# exact-spelling map (kept for compatibility with msp_standarize(LL, MSP_dict) callers)
MSP_dict = {}
for _v, _variants in _MSP_variants.items():
    for _k in _variants:
        for _case in (_k, _k.lower(), _k.upper(), _k.capitalize()):
            if _case != _v:
                MSP_dict.setdefault(_case, _v)


def _norm_key(k):
    """Normalize a field name: lower case, alphanumerics only."""
    return re.sub(r'[^a-z0-9]', '', str(k).lower())

_NORM_VARIANTS = {canon: [_norm_key(v) for v in variants] for canon, variants in _MSP_variants.items()}
_NUM_PEAKS_KEYS = set(_NORM_VARIANTS['Num Peaks'])
_NUMBER = re.compile(r'[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?')


def parse_ri_value(value):
    """
    Retention index value to float, or None.
    Handles plain numbers and forms like 'SemiStdNP=3296.3', '1817.1 (lit)', '1305/2'
    (first number after '=' if present, else first number).
    """
    if value is None:
        return None
    if isinstance(value, (int, float)):
        return None if value != value else float(value)
    s = str(value)
    if '=' in s:
        s = s.split('=', 1)[1]
    m = _NUMBER.search(s)
    return float(m.group()) if m else None


def standardize_entry(e):
    """
    Add canonical keys (see _MSP_variants) to entry dict e, in place; original keys are kept.
    An existing canonical key is not overwritten. 'RI' is converted to float (or None).
    """
    norm_to_key = {}
    for k in e:
        norm_to_key.setdefault(_norm_key(k), k)
    for canon, variants in _NORM_VARIANTS.items():
        if canon in e:
            continue
        for v in variants:
            if v in norm_to_key:
                e[canon] = e[norm_to_key[v]]
                break
    if 'RI' in e:
        e['RI'] = parse_ri_value(e['RI'])
    return e


def msp_standarize(LL, MSP_dict=None):
    """Standardize keys of a list of MSP entry dicts, in place. MSP_dict is ignored
    (kept for backward compatibility); matching uses normalized _MSP_variants."""
    for e in LL:
        standardize_entry(e)
    return LL


def _parse_peak_line(line):
    """
    Peaks from one data line. Handles 'mz int', tabs, quoted annotations
    ('100 45 "loss of H2O"'), several pairs per line separated by ';' or ',' or
    whitespace, '(mz int)' groups and 'mz:int' pairs.
    """
    line = re.sub(r'"[^"]*"', ' ', line)            # drop quoted annotations
    peaks = []
    segments = [s for s in re.split(r'[;()]', line) if s.strip()]
    for seg in segments:
        nums = _NUMBER.findall(seg.replace(':', ' ').replace(',', ' '))
        if len(nums) % 2 == 0:                         # 'mz int' or several pairs on one line
            peaks.extend((float(nums[i]), float(nums[i + 1])) for i in range(0, len(nums), 2))
        elif len(nums) >= 3:                           # 'mz int extra': unquoted 3rd column ignored
            peaks.append((float(nums[0]), float(nums[1])))
    return peaks


def parse_msp_to_listdict(file, field_separator=':', return_peaks=True, standardize=True):
    """
    Parse a msp file within reasonable size.
    Return list of dictionaries, [{key, value pairs}, ...];
    entry['peaks'] contains list of (mz, intensity) tuples if return_peaks.

    Tolerant to common variants:
    - LF or CRLF line endings, blank lines with whitespace, '#' comment lines;
    - records not separated by blank lines (a new 'Name' starts a new record);
    - 'Key: value', 'Key:value', values containing ':' (split at first separator only);
    - repeated keys: first value is kept (e.g. multiple Synon lines);
    - peak lines (starting with a number or '(') in formats handled by _parse_peak_line,
      including lines with annotations that contain ':'.
    With standardize=True (default), canonical keys are added by standardize_entry,
    e.g. 'RI' (float) from 'Retention_index', 'ExactMass' from 'MW', 'Num Peaks' from 'Num.peaks'.
    """
    with open(file, encoding='utf-8', errors='replace') as fh:
        text = fh.read()

    features = []
    d, peaks = {}, []

    def close():
        nonlocal d, peaks
        if d or peaks:
            if return_peaks:
                d['peaks'] = peaks
            if standardize:
                standardize_entry(d)
            features.append(d)
        d, peaks = {}, []

    for raw in text.splitlines():
        line = raw.strip()
        if not line:
            close()
            continue
        if line.startswith('#'):
            continue
        if line[0].isdigit() or line[0] in '.(+-':
            if return_peaks:
                peaks.extend(_parse_peak_line(line))
            continue
        if field_separator in line:
            key, value = line.split(field_separator, 1)
            key, value = key.strip(), value.strip()
            if _norm_key(key) == 'name' and (peaks or any(_norm_key(k) == 'name' for k in d)):
                close()                                # record without blank-line separator
            if key not in d:
                d[key] = value
        elif return_peaks:
            peaks.extend(_parse_peak_line(line))
    close()
    return features


#
# MGF
#

def _parse_mgf_number(value):
    """First number in a string as float, or None."""
    m = _NUMBER.search(str(value))
    return float(m.group()) if m else None


def _parse_mgf_pepmass(value):
    """PEPMASS "mz [intensity [charge]]" -> (mz, intensity); None for missing or non-numeric parts."""
    parts = str(value).replace(',', ' ').split()
    mz = _parse_mgf_number(parts[0]) if parts else None
    intensity = _parse_mgf_number(parts[1]) if len(parts) > 1 else None
    return mz, intensity


def _parse_mgf_charge(value):
    """
    CHARGE to signed int: '2+' -> 2, '2-' -> -2, '-1' -> -1, '+3' -> 3, '1' -> 1.
    For multiple charges ('2+ and 3+', '2+,3+') the first one is used. None if unparsable.
    """
    m = re.search(r'([+-]?)\s*(\d+)\s*([+-]?)', str(value))
    if not m:
        return None
    sign = -1 if '-' in (m.group(1) + m.group(3)) else 1
    return sign * int(m.group(2))


def parse_mgf(file_path, standardize=True, warn=True):
    """
    Parse an MGF file. Returns a list of spectrum dicts:
        'params'      : header fields, keys upper-cased; PEPMASS as (mz, intensity),
                        CHARGE as signed int (e.g. '2-' -> -2)
        'mzs', 'intensities' : peak lists (as before)
        'peaks'       : list of (mz, intensity) tuples, consistent with parse_msp_to_listdict
        'precursor_mz', 'precursor_intensity', 'charge'
        'rtime'       : retention time in seconds, from RTINSECONDS or RTINMINUTES (else None)
    With standardize=True, canonical keys (Name, InChIKey, SMILES, RI, ...; see _MSP_variants)
    are also added at top level from the params.

    Tolerant to: LF/CRLF, lower-case 'begin ions', comment lines starting with
    '#', ';', '!' or '/', blank lines, values containing '=', peaks separated by
    whitespace, tabs or comma, annotation columns after intensity, non-numeric PEPMASS.
    A spectrum missing 'END IONS' (followed by a new 'BEGIN IONS' or end of file)
    is kept, with a warning if warn=True. Malformed peak lines are skipped.
    """
    spectra = []
    current = None
    n_unterminated = 0

    def finish(spec):
        params = spec['params']
        pep = params.get('PEPMASS')
        spec['precursor_mz'], spec['precursor_intensity'] = pep if pep else (None, None)
        spec['charge'] = params.get('CHARGE')
        rtime = None
        if params.get('RTINSECONDS') not in (None, ''):
            rtime = _parse_mgf_number(params['RTINSECONDS'])
        elif params.get('RTINMINUTES') not in (None, ''):
            rt_min = _parse_mgf_number(params['RTINMINUTES'])
            rtime = None if rt_min is None else rt_min * 60
        spec['rtime'] = rtime
        spec['peaks'] = list(zip(spec['mzs'], spec['intensities']))
        if standardize:
            canon = standardize_entry(dict(params))
            for k in _MSP_variants:
                if k in canon and k not in spec and k not in ('Num Peaks', 'PrecursorMZ'):
                    spec[k] = canon[k]
        spectra.append(spec)

    with open(file_path, encoding='utf-8', errors='replace') as f:
        for raw in f:
            line = raw.strip()
            if not line or line[0] in '#;!/':
                continue
            upper = line.upper()

            if upper == 'BEGIN IONS':
                if current is not None:          # previous spectrum without END IONS
                    n_unterminated += 1
                    finish(current)
                current = {'params': {}, 'mzs': [], 'intensities': []}

            elif upper == 'END IONS':
                if current is not None:
                    finish(current)
                    current = None

            elif current is not None:
                if line[0].isdigit() or line[0] in '.+-':
                    nums = _NUMBER.findall(line.replace(',', ' ').split('"')[0])
                    if len(nums) >= 2:
                        current['mzs'].append(float(nums[0]))
                        current['intensities'].append(float(nums[1]))
                elif '=' in line:
                    key, value = line.split('=', 1)
                    key, value = key.strip().upper(), value.strip()
                    if key == 'PEPMASS':
                        value = _parse_mgf_pepmass(value)
                    elif key == 'CHARGE':
                        value = _parse_mgf_charge(value)
                    current['params'][key] = value

    if current is not None:                      # last spectrum without END IONS
        n_unterminated += 1
        finish(current)
    if warn and n_unterminated:
        print(f"parse_mgf: {n_unterminated} spectra without 'END IONS' in {file_path}; kept.")
    return spectra
