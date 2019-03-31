#### NOTE: Before calling these functions,
# in your main Python file, you must:
# cimport numpy as np
# np.import_array()


from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool
from libc.stdint cimport uint8_t
from libc.string cimport memcpy

from cpython cimport PyObject,Py_DECREF
cimport numpy as np
from numpy cimport npy_intp, PyArray_SimpleNewFromData, PyArray_SimpleNew

cdef extern from "ncorr.h" namespace "ncorr::Data2D":
    ctypedef ptrdiff_t difference_type

cdef extern from "ncorr.h" namespace "ncorr" nogil:

    cdef cppclass Array2D[T]:
        Array2D()
        Array2D(size_t h,size_t w,T val)
        T& operator()(difference_type,difference_type)

        difference_type width()
        difference_type height()

        T *get_pointer()
        
        Array2D[bool] operator>(double) #Needed for the input roi image >
        pass

    cdef cppclass ROI2D: #Needed to convert Array2D[bool] to python object
        ROI2D()
        ROI2D(Array2D[bool], difference_type diff) #Leads to an indent level error
        const Array2D[bool]& get_mask()
        
        pass

    cdef cppclass Data2D:
        difference_type data_height()
        difference_type data_width()
        Array2D[double]& get_array() 
        const ROI2D& get_roi() 
        difference_type get_scalefactor()


        pass

    cdef cppclass Disp2D:
        difference_type data_height()
        difference_type data_width()
        Data2D& get_v()
        Data2D& get_u()
        ROI2D& get_roi()
        difference_type get_scalefactor()
        pass


    cdef cppclass DIC_analysis_output:
        vector[Disp2D] disps

        pass
        
    cdef cppclass strain_analysis_output:
        pass
    cdef cppclass strain_analysis_input:
        pass
    cdef cppclass Image2D:
        Image2D(string) #Required to prevent crash in optimizebuiltincalls, constructor in C++, Converting Python object to 'Array2d[bool]'
        Image2D(Array2D[double] &)
        Array2D[double] get_gs() #Needed to object 'Image2D', converting python object
        pass
    cdef enum INTERP: #All of this is needed for undeclared name not builtin
        NEAREST "ncorr::INTERP::NEAREST"
        LINEAR "ncorr::INTERP::LINEAR"
        QUINTIC_BSPLINE_PRECOMPUTE "ncorr::INTERP::QUINTIC_BSPLINE_PRECOMPUTE"
        # ...
        pass
    cdef enum SUBREGION:
        CIRCLE "ncorr::SUBREGION::CIRCLE"
        pass
    cdef enum DIC_analysis_config:
        NO_UPDATE "ncorr::DIC_analysis_config::NO_UPDATE"
        KEEP_MOST_POINTS "ncorr::DIC_analysis_config::KEEP_MOST_POINTS"
        REMOVE_BAD_POINTS "ncorr::DIC_analysis_config::REMOVE_BAD_POINTS"
        pass
        
    
    
    cdef cppclass DIC_analysis_input:
        DIC_analysis_input()
        DIC_analysis_input(const vector[Image2D]&, const ROI2D&, difference_type, INTERP, SUBREGION, difference_type, difference_type, DIC_analysis_config, bool)
        @staticmethod
        DIC_analysis_input load(string name)
        
        
        
        pass
    cdef void Array2Dbool_save "save" (const Array2D[bool] &A, const string &filename)
    cdef void Array2Ddouble_save "save" (const Array2D[double] &A, const string &filename)
        
    cdef void DIC_analysis_output_save "save" (const DIC_analysis_output&, const string&)	
    
    cdef DIC_analysis_output DIC_analysis(const DIC_analysis_input&)
        
    cdef DIC_analysis_output set_units(const DIC_analysis_output&, const string&, double)
    
    
    #cdef save(const DIC_analysis_input&, const string&)
    #pass
#    
#    cdef save(const DIC_analysis_output&, const string&);
#    pass
    
    
    pass
    

cdef inline doublearray2numpy(const Array2D[double] &array):
    cdef npy_intp[2] dims
    cdef double *arraydata
    #cdef PyObject *newarray_c
    dims[0]=array.width()
    dims[1]=array.height()
    # Array2D is indexed, (rows, columns), Fortran-style (y,x)
    # ... we ready it C-style (x,y) then transpose

    
    arraydata = array.get_pointer()
    
    newarray = PyArray_SimpleNewFromData(2,dims,np.NPY_DOUBLE,arraydata).T
    
    #newarray=<object>newarray_c
    #Py_DECREF(newarray_c)
    
    return newarray.copy() # keep our own copy so the data is properly reference counted

cdef inline boolarray2numpy(const Array2D[bool] &array):
    cdef npy_intp[2] dims
    cdef bool *arraydata
    cdef size_t nrows
    cdef size_t ncols
    cdef size_t rowcnt
    cdef size_t colcnt
    cdef uint8_t *dataptr
    cdef np.ndarray[np.uint8_t,ndim=2,mode="fortran"] newarray
    
    #cdef PyObject *newarray_c
    dims[0]=array.width()
    dims[1]=array.height()
    # Array2D is indexed, (rows, columns), Fortran-style (y,x)
    # ... we create it C-style (x,y) then transpose

    nrows=array.height()
    ncols=array.width()

    
    arraydata = array.get_pointer()
    
    newarray = PyArray_SimpleNew(2,dims,np.NPY_UINT8).T
    dataptr=<uint8_t *>newarray.data
    
    with nogil: 
       for rowcnt in range(nrows):
           for colcnt in range(ncols):
               dataptr[colcnt*nrows + rowcnt] = arraydata[colcnt*nrows + rowcnt]
               pass
           pass
       pass
    #newarray=<object>newarray_c
    #Py_DECREF(newarray_c)
    
    return newarray


cdef inline Array2D[double] numpy2doublearray(np.ndarray[np.float64_t,ndim=2,mode="fortran"] nparray):
    cdef Array2D[double] resarray
    cdef double *res_ptr
    # Array2D is indexed, (rows, columns), Fortran-style (y,x)
    resarray=Array2D[double](nparray.shape[0],nparray.shape[1],0.0)
    res_ptr = resarray.get_pointer()

    memcpy(res_ptr,nparray.data,sizeof(double)*nparray.shape[0]*nparray.shape[1])
    
    
    return resarray

cdef inline Array2D[bool] numpy2boolarray(np.ndarray[np.uint8_t,ndim=2,mode="fortran"] nparray):
    cdef Array2D[bool] resarray
    cdef bool *res_ptr
    cdef size_t nrows
    cdef size_t ncols

    cdef size_t rowcnt
    cdef size_t colcnt
    cdef uint8_t *dataptr
    
    nrows=nparray.shape[0]
    ncols=nparray.shape[1]

    dataptr=<uint8_t *>nparray.data
    
    # Array2D is indexed, (rows, columns), Fortran-style (y,x)
    resarray=Array2D[bool](nrows,ncols,False)
    res_ptr = resarray.get_pointer()

    with nogil: 
       for rowcnt in range(nrows):
           for colcnt in range(ncols):
               res_ptr[colcnt*nrows + rowcnt] = dataptr[colcnt*nrows + rowcnt]
               pass
           pass
       pass
    #Array2Dbool_save(resarray,"/tmp/array2dbool.txt")
    
    return resarray

