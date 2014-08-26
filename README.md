# VICUS

VICUS is a [linear programming](https://en.wikipedia.org/wiki/Linear_programming) optimisation model for capacity expansion planning and unit commitment for localised energy systems. Its name, latin for village, is derived from its bigger sister model [URBS](https://github.com/tum-ens/urbs) (latin for city). Its name stems from its origin as a model for optimisation for urban energy systems.

## Features

  * VICUS is a linear programming model for multi-commodity energy systems with a focus on optimal storage sizing and use.
  * It finds the minimum cost energy system to satisfy given demand timeseries for possibly multiple commodities (e.g. electricity).
  * By default, operates on hourly-spaced timesteps (configurable).
  * Thanks to [pandas](https://pandas.pydata.org), data analysis code is short and clean.
  * The model itself is quite small (<40 kB source code) thanks to relying on the [Coopr](https://software.sandia.gov/trac/coopr)/[Pyomo](https://software.sandia.gov/trac/coopr/wiki/Pyomo) and includes reporting and plotting functionality.

## Installation

### Windows

For all packages, best take the latest release or release candidate version. Both 32 bit and 64 bit versions work, though 64 bit is recommended.

  1. **[Python 2.7](https://python.org/download)**. Python 3 support is not possible yet, but planned once all used packages support it.
  2. **[pip](https://pip.pypa.io/en/latest/installing.html)**.The Python package manager. It allows to install many Python packages with a simple command. 
      1. After installation, add `C:\Python27\Scripts` to environment variable "Path" ([how](http://geekswithblogs.net/renso/archive/2009/10/21/how-to-set-the-windows-path-in-windows-7.aspx)), so that the `pip` command becomes available on the command prompt.
  3. **IPython**: execute `pip install ipython` in a command prompt.
  4. **SciPy stack:** These require binary installers, made available and maintained by [C. Gohlke](http://www.lfd.uci.edu/~gohlke/pythonlibs/). *How to select the correct file:* Download the newest stable version of each package, whose filename suffix matches both "bitness" (32 bit or 64 bit) and Python version (i.e. 2.7).  
      1. [NumPy](http://www.lfd.uci.edu/~gohlke/pythonlibs/#numpy)
      2. [SciPy](http://www.lfd.uci.edu/~gohlke/pythonlibs/#scipy)
      3. [matplotlib](http://www.lfd.uci.edu/~gohlke/pythonlibs/#matplotlib), requires [dateutil](http://www.lfd.uci.edu/~gohlke/pythonlibs/#python-dateutil), [pytz](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pytz), [pyparsing](http://www.lfd.uci.edu/~gohlke/pythonlibs/#pyparsing) and [six](http://www.lfd.uci.edu/~gohlke/pythonlibs/#six). 
      4. As a test, you can try start `ipython --pylab` and have a MATLAB-style command line with plotting capabilities. If you receive message about "ipython could not be found", check if the `C:\Python27\Scripts` is added to the "Path" system variable as described in step 2.i. above.
  5. **[pandas](https://pypi.python.org/pypi/pandas#downloads)**: its [Series](http://pandas.pydata.org/pandas-docs/stable/dsintro.html#series) and [DataFrame](http://pandas.pydata.org/pandas-docs/stable/dsintro.html#dataframe) are used for representing all model input and output. Its capabilities are exploited to write short analysis scripts in `runme.py` and `comp.py`, as well as in the functions `vicus.plot` and `vicus.report`.
  6. **[Coopr](https://software.sandia.gov/trac/coopr/downloader/)**: minimum version 3.5 or the VOTD (Version of the Day) installer. As of 2014-08-01, only the latter is available for Windows users.
  7. **Solver**: [GLPK](http://winglpk.sourceforge.net/). 
      1. Simply unzip the latest version somewhere, e.g. `C:\GLPK`. 
      2. Then add the subdirectory `w64`, which contains `glpsol.exe`, to the system path (like in step 2.i.), so that the `glpsol` command is available on the command prompt.
  8. **Excel** reading/writing: `pip install xlrd xlwt openpyxl==1.8.6`

### Linux

Use the package manager to get all the packages listed in the Windows installation section. Below is the installation procedure for Ubuntu & Debian. Other distributions might have slightly different package names or differing procedures to get the individual packages to run:

  - **Everything** except Coopr & Excel I/O `sudo apt-get install python python-pip python-numpy python-scipy python-matplotlib ipython ipython-notebook python-pandas python-sympy python-nose glpk-utils`
  - **Coopr & Excel I/O** `sudo pip install coopr xlwt xlrd openpyxl==1.8.6`

## Get started

Once installation is complete, clone (or download) this repository and execute the runme script on the command prompt:

    git clone https://github.com/ojdo/vicus.git
    cd vicus
    python runme.py

About a minute later, two pictures `plot-ElecAC.png`, `plot-ElecDC.png` and a spreadsheet `report.xlsx` should have been generated.

## Next steps

  1. Read the source code of `runme.py`.
  2. Quickly scan through `vicus.py`, read function names and their docstrings.
  3. Fire up IPython (`ipython --pylab`) and run the script from there using the run command: `run runme`. Then use `whos` and inspect the workspace afterwards. See what you can do (analyses, plotting) with the DataFrames returned by `vicus.get_constants` and `vicus.get_timeseries` after a successful optimisation run. Take `vicus.plot` and `vicus.report` functions as inspriation and the [pandas docs](http://pandas.pydata.org/pandas-docs/stable/) as reference.
  
## Further reading

  - The book [Python for Data Analysis](http://shop.oreilly.com/product/0636920023784.do) best summarises the capabilities of the packages installed here. It starts with IPython, then adds NumPy, slowly fades to pandas and then shows first basic, then advanced data conversion and analysis recipes. Visualisation with matplotlib is given its own chapter, both with and without pandas.
  - For a huge buffet of appetizers, showing the capabilities of the Python for scientific computing, I recommend browsing this [gallery of interesting IPython Notebooks](https://github.com/ipython/ipython/wiki/A-gallery-of-interesting-IPython-Notebooks).
