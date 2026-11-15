"""
Placeholder for m/z calibration of features against a known compound database.

Moved out of LCMS_Annotation (annotate.py) in asari 2.0, not used in processing or annotation.
Example use:

    from jms.dbStructures import knownCompoundDatabase
    from asari.tools.db_mass_calibrate import db_mass_calibrate

    KCD = knownCompoundDatabase()       # with mass_indexed_compounds and emp_cpds_trees loaded
    list_features.sort(key=lambda x: x['peak_area'], reverse=True)
    db_mass_calibrate(list_features, KCD, mode='pos')
"""


def db_mass_calibrate(list_features, KCD, mode='pos',
                      max_features=1000, required_calibrate_threshold=0.000002, mz_tolerance_ppm=10):
    """
    Use KCD.evaluate_mass_accuracy_ratio to check systematic mass shift,
    which is calculated as the average ppm difference between measured m/z and theoretical values.
    If greater than required_calibrate_threshold (default 2 ppm),
    calibrate m/z values of all features in list_features (in place).

    Parameters
    ----------
    list_features : list
        features as dicts with 'mz', e.g. from read_features_from_asari_table,
        sorted by intensity (decreasing) so that the top max_features are used as landmarks.
    KCD : jms.dbStructures.knownCompoundDatabase
        known compound database with mass_indexed_compounds loaded.
    mode : str, optional, default: 'pos'
        ionization mode.
    max_features : int, optional, default: 1000
        number of top features used to estimate mass shift.
    required_calibrate_threshold : float, optional, default: 0.000002
        if the mass shift exceeds this value, mass correction will be applied.
    mz_tolerance_ppm : float, optional, default: 10
        m/z tolerance to match features to database.

    Returns
    -------
    mass_accuracy_ratio, or None if too few features matched to the database.
    Corrected features have 'mz' updated and 'mz_corrected_by_division' added.
    """
    mz_landmarks = [f['mz'] for f in list_features[:max_features]]
    mass_accuracy_ratio = KCD.evaluate_mass_accuracy_ratio(
        mz_landmarks, mode=mode, mz_tolerance_ppm=mz_tolerance_ppm)
    if mass_accuracy_ratio:
        if abs(mass_accuracy_ratio) > required_calibrate_threshold:
            print("Mass shift is greater than %2.1f ppm. Correction applied."
                  %(required_calibrate_threshold*1000000))
            _correction = mass_accuracy_ratio + 1
            for F in list_features:
                F['mz'] = F['mz'] / _correction
                F['mz_corrected_by_division'] = _correction
    else:
        print("Mass accuracy check is skipped, too few of %d features matched to database."
              %len(mz_landmarks))
    return mass_accuracy_ratio
