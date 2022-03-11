### to add to doc/

notebook explaining the algorithms.

notebook example of how to use the library for advanced functions.

CenturionTree
=============
a dictionary, indexing mzList by 100*mz bins.
Because most high-resolution mass spectrometers measure well under 0.01 amu, 
one only needs to search the corresponding 0.01 bin and two adjacent bins (to capture bordering values).

Use scan numbers 
================
Because scan numbers are integers, they are efficient as indices and should be used for most low-level operations.
When real retention time is used, they are float numbers and not good for indices, 
requiring many more comparison operations and decreasing performance.

    
To improve performance in next version:

Use C to rewrite chromatogram constructor.
After initial samples, the peak detection of most features can start from building EIC in C, 
to reduce workload in scipy.signal.find_peak.



Notebooks 
=========

- Single sample processing, inspection, and determine ppm precision.

- Process data without upfront LC alignment

- Annotate and search

- OpenMS based workflow


- batched processing; separate CMAP per processing batch and assemble afterwards. 
    Or, how to combine multiple asari projects (based on empCpds and massTracks).



