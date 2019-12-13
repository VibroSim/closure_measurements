import shutil
import re
import os
import numpy as np
import os.path
import subprocess
from setuptools import setup
from setuptools.command.install_lib import install_lib
from setuptools.command.install import install
from setuptools.command.build_ext import build_ext
import setuptools.command.bdist_egg
import sys
#import glob

from Cython.Build import cythonize


extra_compile_args = {
    "msvc": ["/std:c++14"],
    #"unix": ["-O0", "-g", "-Wno-uninitialized"),    # Replace the line below with this line to enable debugging of the compiled extension
    "unix": ["-O5","-Wno-uninitialized","-std=c++11"],
    #"unix": ["-O0","-Wno-uninitialized","-std=c++11","-g"],
    "clang": ["-O5","-Wno-uninitialized"],
}

extra_include_dirs = {
    "msvc": ["c:\\suitesparse", np.get_include() ],
    "unix": ["/usr/include/suitesparse",np.get_include() ],
    "clang": ["/usr/include/suitesparse", np.get_include() ],
}

extra_libraries = {
    "msvc": ["ncorr","opencv_core","opencv_imgproc","opencv_videoio","opecv_highgui","cholmod","spqr","fftw3"],
    "unix": ["ncorr","opencv_core","opencv_imgproc","opencv_videoio","opecv_highgui","cholmod","spqr","fftw3"],
    "clang": ["ncorr","opencv_core","opencv_imgproc","opencv_videoio","opecv_highgui","cholmod","spqr","fftw3"],
}

extra_link_args = {
    "msvc": [],
    "unix": [],
    "clang": [],
}

class build_ext_compile_args(build_ext):
    def build_extensions(self):
        compiler=self.compiler.compiler_type
        for ext in self.extensions:
            if compiler in extra_compile_args:
                ext.extra_compile_args=extra_compile_args[compiler]
                ext.extra_link_args=extra_link_args[compiler]
                ext.include_dirs.extend(list(extra_include_dirs[compiler]))
                ext.libraries.extend(list(extra_libraries[compiler]))
                pass
            else:
                # use unix parameters as default
                ext.extra_compile_args=extra_compile_args["unix"]
                ext.extra_link_args=extra_link_args["unix"]
                ext.include_dirs.extend(list(extra_include_dirs["unix"]))
                ext.libraries.extend(extra_libraries["unix"])
                pass
                
            pass
            
        
        build_ext.build_extensions(self)
        pass
    pass


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
        versionraw=subprocess.check_output(["git","describe","--tags","--match=v*"],stderr=subprocess.STDOUT).decode('utf-8').strip()
        # versionraw is like v0.1.0-50-g434343
        # for compatibility with PEP 440, change it to
        # something like 0.1.0+50.g434343
        matchobj=re.match(r"""v([^.]+[.][^.]+[.][^-.]+)(-.*)?""",versionraw)
        version=matchobj.group(1)
        if matchobj.group(2) is not None:
            version += '+'+matchobj.group(2)[1:].replace("-",".")
            pass
        pass
    except subprocess.CalledProcessError:
        # Ignore error, falling back to above version string
        pass

    if modified and version.find('+') >= 0:
        version += ".modified"
        pass
    elif modified:
        version += "+modified"
        pass
    pass
else:
    version = "UNKNOWN"
    pass

print("version = %s" % (version))

closure_measurements_package_files=[ "qagse_fparams.c","pt_steps/*"]


console_scripts=["closure_measurement_dic","closure_measurement_coords","closure_measurement_processing"]
#gui_scripts = []  # Could move graphical scripts into here to eliminate stdio window on Windows (where would error messages go?)

console_scripts_entrypoints = [ "%s = closure_measurements.bin.%s:main" % (script,script.replace("-","_")) for script in console_scripts ]

#gui_scripts_entrypoints = [ "%s = limatix.bin.%s:main" % (script,script.replace("-","_")) for script in gui_scripts ]


ext_modules=cythonize("closure_measurements/*.pyx",language="c++")

emdict=dict([ (module.name,module) for module in ext_modules])

#correlate_pyx_ext=emdict['closure_measurements.correlate']
##correlate_pyx_ext.sources.append("extra_c_file.c")
#correlate_pyx_ext.extra_compile_args=['-std=c++11', '-O0','-g']
#correlate_pyx_ext.include_dirs = ['/usr/include/suitesparse']
#correlate_pyx_ext.libraries = ['ncorr','opencv_core','opencv_imgproc','opencv_videoio','opencv_highgui','cholmod','spqr','fftw3']

##correlate_pyx_ext.extra_compile_args=['-fopenmp','-O3']
##correlate_pyx_ext.extra_link_args=['-lgomp']



setup(name="closure_measurements",
      description="Inversion of crack closure from DIC measurements",
      author="Chris Giuffre and Stephen D. Holland",
      version=version,
      # url="http://",
      zip_safe=False,
      packages=["closure_measurements",
                "closure_measurements.bin"],
      cmdclass={"install_lib": install_lib_save_version,
                "build_ext": build_ext_compile_args},
      #data_files=[ ("share/closure_measurements/pt_steps",pt_steps_files),]
      package_data={"closure_measurements": closure_measurements_package_files},
      entry_points={ "limatix.processtrak.step_url_search_path": [ "limatix.share.pt_steps = closure_measurements:getstepurlpath" ],
                    "console_scripts": console_scripts_entrypoints,
                    #"gui_scripts": gui_scripts_entrypoints          
                },
      ext_modules=ext_modules)


