# VICUS

VICUS is a [linear programming](https://en.wikipedia.org/wiki/Linear_programming) optimisation model for capacity expansion planning and unit commitment for localised energy systems. Its name, latin for village, is derived from its bigger sister model [URBS](https://github.com/tum-ens/urbs) (latin for city). Its name stems from its origin as a model for optimisation for urban energy systems.

## Features

  * VICUS is a linear programming model for multi-commodity energy systems with a focus on optimal storage sizing and use.
  * It finds the minimum cost energy system to satisfy given demand timeseries for possibly multiple commodities (e.g. electricity).
  * By default, operates on hourly-spaced timesteps (configurable).
  * Thanks to [pandas](https://pandas.pydata.org), data analysis code is short and clean.
  * The model itself is quite small (<40 kB source code) thanks to relying on the [Coopr](https://software.sandia.gov/trac/coopr)/[Pyomo](https://software.sandia.gov/trac/coopr/wiki/Pyomo) and includes reporting and plotting functionality.

## Installation

Follow the installation instructions for the [urbs](https://github.com/tum-ens/urbs#installation), then continue with the section **Get started** below.

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
