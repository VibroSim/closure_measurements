import sys
import os
import os.path
import copy
import ast
import posixpath
import numpy as np
import subprocess
import limatix.timestamp




from limatix import dc_value
from limatix.dc_value import hrefvalue as hrefv    
from limatix import xmldoc

if sys.version_info[0] < 3:
    input = raw_input   # backwards compatibility with Python 2.x
    pass

def run(_xmldoc,_element,dc_closureprofileplot_href):

    tipfit_side1 = _xmldoc.xpathsinglecontext(_element,"dc:dic_tip_fit[@side='1']",default=None)
    tipfit_side2 = _xmldoc.xpathsinglecontext(_element,"dc:dic_tip_fit[@side='2']",default=None)

    closureplot_side1 = _xmldoc.xpathsinglecontext(_element,"dc:closureplot_side1",default=None)
    closureplot_side2 = _xmldoc.xpathsinglecontext(_element,"dc:closureplot_side2",default=None)

    first_href = None

    print("View the following images: ")

    

    if tipfit_side1 is not None:
        tipfit_side1_href = hrefv.fromxml(_xmldoc,tipfit_side1)
        print("  %s" % (os.path.split(tipfit_side1_href.getpath())[1]))
        first_href = tipfit_side1_href
        pass

    if tipfit_side2 is not None:
        tipfit_side2_href = hrefv.fromxml(_xmldoc,tipfit_side2)
        print("  %s" % (os.path.split(tipfit_side2_href.getpath())[1]))

        if first_href is None:
            first_href = tipfit_side2_href
            pass
        pass

    print("  %s" % (os.path.split(dc_closureprofileplot_href.getpath())[1]))

    
    if closureplot_side1 is not None:
        closureplot_side1_href = hrefv.fromxml(_xmldoc,closureplot_side1)
        print("  %s" % (os.path.split(closureplot_side1_href.getpath())[1]))
        pass

    if closureplot_side2 is not None:
        closureplot_side2_href = hrefv.fromxml(_xmldoc,closureplot_side2)
        print("  %s" % (os.path.split(closureplot_side2_href.getpath())[1]))
        pass
        
        

    if first_href is None:
        raise ValueError("No measurements found!")
        
    viewerproc = subprocess.Popen(["geeqie",first_href.getpath()])

    print("Verify:")
    print("  * That tip fits are reasonable")
    print("  * That the resulting closure profile is reasonable")
    print("  * That the resulting closure plots are reasonable")

    print("If not, answer 'n' and make appropriate corrections if possible.")


    OK_text = input("All are OK [y/N]? ")      
    if OK_text.lower()=='y' or OK_text.lower()=="yes":
        ret = {"dc:ClosureQA": limatix.timestamp.now().isoformat()}
        pass
    else:
        raise ValueError("Closure QA failed")
        
    return ret
            
