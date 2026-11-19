# requires pyopenms, not in Asari dependency

import sys, os, glob
from multiprocessing import Pool
import pyopenms as oms

src, dst = sys.argv[1], sys.argv[2]

def conv(f):
    out = os.path.join(dst, os.path.basename(f)[:-6] + '.mzML')
    if os.path.exists(out):
        return f, 'exists'
    exp = oms.MSExperiment()
    oms.MzXMLFile().load(f, exp)
    n_profile = sum(1 for s in exp if s.getType() == oms.SpectrumSettings.SpectrumType.PROFILE)
    if n_profile:
        oms.PeakPickerHiRes().pickExperiment(exp, exp_out := oms.MSExperiment(), True)
        exp = exp_out
    oms.MzMLFile().store(out + '.tmp', exp)
    os.rename(out + '.tmp', out)
    return f, f'{exp.getNrSpectra()} spectra, picked={bool(n_profile)}'

if __name__ == '__main__':
    files = sorted(glob.glob(os.path.join(src, '*.mzXML')))
    if len(sys.argv) > 3: files = files[:int(sys.argv[3])]
    with Pool(6) as p:
        for f, msg in p.imap_unordered(conv, files):
            print(os.path.basename(f), msg, flush=True)
