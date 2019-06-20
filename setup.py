import shutil
import os
import os.path
import subprocess
from setuptools import setup
from setuptools.command.install_lib import install_lib
from setuptools.command.install import install
import setuptools.command.bdist_egg
import sys
#import glob

from Cython.Build import cythonize

class install_lib_save_version(install_lib):
    """Save version information"""
    def run(self):
        install_lib.run(self)
        
        for package in self.distribution.command_obj["build_py"].packages:
            install_dir=os.path.join(*([self.install_dir] + package.split('.')))
            fh=open(os.path.join(install_dir,"version.txt"),"w")
            fh.write("%s\n" % (version))  # version global, as created below
            fh.close()
            pass
        pass
    pass


# Extract GIT version
if os.path.exists(".git"):
    # Check if tree has been modified
    modified = subprocess.call(["git","diff-index","--quiet","HEAD","--"]) != 0
    
    gitrev = subprocess.check_output(["git","rev-parse","HEAD"]).strip()

    version = "git-%s" % (gitrev)

    # See if we can get a more meaningful description from "git describe"
    try:
        version=subprocess.check_output(["git","describe"]).strip()
        pass
    except subprocess.CalledProcessError:
        # Ignore error, falling back to above version string
        pass

    if modified:
        version += "-MODIFIED"
        pass
    pass
else:
    version = "UNKNOWN"
    pass

closure_measurements_package_files=[ "qagse_fparams.c","pt_steps/*"]


console_scripts=["closure_measurement_dic","closure_measurement_coords","closure_measurement_processing"]
#gui_scripts = []  # Could move graphical scripts into here to eliminate stdio window on Windows (where would error messages go?)

console_scripts_entrypoints = [ "%s = closure_measurements.bin.%s:main" % (script,script.replace("-","_")) for script in console_scripts ]

#gui_scripts_entrypoints = [ "%s = limatix.bin.%s:main" % (script,script.replace("-","_")) for script in gui_scripts ]


ext_modules=cythonize("closure_measurements/*.pyx",language="c++")

emdict=dict([ (module.name,module) for module in ext_modules])

correlate_pyx_ext=emdict['closure_measurements.correlate']
#correlate_pyx_ext.sources.append("extra_c_file.c")
correlate_pyx_ext.extra_compile_args=['-std=c++11', '-O0','-g']
correlate_pyx_ext.include_dirs = ['/usr/include/suitesparse']
correlate_pyx_ext.libraries = ['ncorr','opencv_core','opencv_imgproc','opencv_videoio','opencv_highgui','cholmod','spqr','fftw3']

#correlate_pyx_ext.extra_compile_args=['-fopenmp','-O3']
#correlate_pyx_ext.extra_link_args=['-lgomp']



setup(name="closure_measurements",
      description="Inversion of crack closure from DIC measurements",
      author="Chris Giuffre and Stephen D. Holland",
      # url="http://",
      zip_safe=False,
      packages=["closure_measurements",
                "closure_measurements.bin"],
      #data_files=[ ("share/closure_measurements/pt_steps",pt_steps_files),]
      package_data={"closure_measurements": closure_measurements_package_files},
      entry_points={ "limatix.processtrak.step_url_search_path": [ "limatix.share.pt_steps = closure_measurements:getstepurlpath" ],
                    "console_scripts": console_scripts_entrypoints,
                    #"gui_scripts": gui_scripts_entrypoints          
                },
      ext_modules=ext_modules)


