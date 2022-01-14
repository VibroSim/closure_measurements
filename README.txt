closure_measurements
--------------------

This Python package attempts to intuit crack closure states from
measured digital image correlation data. 

Requirements
------------
(older versions may work but have not been tested)
  * Python 2.7.5 or 3.4.9 or newer. 
  * scipy 1.0.0 or newer
  * numpy 1.14.3 or newer
  * matplotlib 1.5.3 or newer
  * IPython/Jupyter (recommended)
  * git 2.17.1 or newer 
    (you may delete the .git directory if you prefer not to use version
    control.)
  * crackclosuresim2 v0.4.1 or newer
  * crackclosuresim2
  * Cython and a C compiler. See the README.txt
    with crackclosuresim2 for more information. 
  * Forked version of ncorr from https://github.com/sdh4/ncorr_2D_cpp

On Linux these components, with the exception of crackclosuresim2
and ncorr are usually available as packages from your operating system vendor. 

On Windows/Macintosh it is usually easiest to use a Python distribution 
such as Anaconda https://www.anaconda.com or Canopy 
https://www.enthought.com/product/canopy/ 

These distributions typically provide the full 
Python/Numpy/Matplotlib/IPython stack by default, so you only need
a few more pieces such as Cython, git, and the C compiler. 
64-bit versions of the distributions are recommended

Installing closure_measurements
--------------------------------
From a terminal, command prompt, Anaconda or Canopy terminal, etc. 
change to the angled_friction_model source directory and type:
  python setup.py build
  python setup.py install

Depending on your configuration the 'install' step might require
root or administrator permissions. You can also explicitly specify 
a different Python version or binary. 

Running angled_friction_model
-----------------------------

Try the examples in the 'demos/' subdirectory. 
   e.g. python simple_afm_demo.py

We recommend using an IPython/Jupyter 
Qt console or similar. Usually you will want to 
start your session by initializing matplotlib mode: 
  %matplotlib qt

Then run one of the demos:
  %run test_dic
  %run test_process_dic

When writing your own Python code, you can import the closure_measurements package
with: 
  import closure_measurements


Sponsor Acknowledgment
----------------------

This material is based on work supported by the Air Force Research
Laboratory under Contract #FA8650-16-C-5238 and performed at Iowa
State University

This material is based upon work supported by the Department of the
Air Force under Awards No. FD20301933322 and FD20301933313, Air Force
SBIRs to Core Parts, LLC.

AFRL Public Release Case #AFRL-2021-3480

Any opinions, findings, and conclusions or recommendations expressed
in this material are those of the author(s) and do not necessarily
reflect views of the Department of the Air Force or Core Parts, LLC.

