# IAAFT functions

This repository provides **Matlab** and **Python** functions implementing the IAAFT algorithm for generation of surrogate time series.

The **IAAFT** (**i**terative **a**mplitude **a**djusted **F**ourier **t**ransform) algorithm was first proposed by T. Schreiber and A. Schmitz in Phys. Rev. Lett. **77** (1996), 635.

## License

[![License: LGPL v3](https://img.shields.io/badge/License-LGPL%20v3-blue.svg?style=flat-square)](/LICENSE)

This package is free software. It is distributed under the terms of the GNU General Public License (GPL), version 3.0 - see the LICENSE.txt file for details.

## Authors

- Alessio Perinelli (1), alessio.perinelli@unitn.it
- Leonardo Ricci (1,2), leonardo.ricci@unitn.it

(1) Department of Physics, University of Trento, Trento, Italy  
(2) CIMeC, Center for Mind/Brain Sciences, University of Trento, Rovereto, Italy  

## Matlab
The Matlab implementation of the function was developed and tested on Matlab version R2019a, but should be compatible with any reasonably recent version of Matlab. The package does not require any setup. To make the function available in Matlab, the source (`.m`) file has to be copied in a directory that is part of Matlab _path_. For further details, please refer to the [official Mathworks website/user guide](https://www.mathworks.com/help/matlab/matlab_env/what-is-the-matlab-search-path.html).

## Python
The Python implementation requires Python 3.5 and NumPy 1.15 or later. The package might not work properly with older versions of Python or NumPy.
