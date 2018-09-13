import numpy as np
cimport numpy as np

## References the functions defined in the header of the C library
cdef extern from "multiply.h":
    void multiply_by_10_in_C(double arr[], unsigned int n)


## Functions alone can be used for operations that input
## and output only Python-compatible datatypes
def py_multiply_by_10_in_C(arr) -> None:
    if not arr.flags['C_CONTIGUOUS']:
        ## Makes a contiguous copy of the numpy array
        arr = np.ascontiguousarray(arr) 

    ## Should better explain what's happening here - takes a numpy 'view'
    ## of the array I think, using special syntax
    cdef double[::1] arr_memview = arr

    multiply_by_10_in_C(&arr_memview[0], arr_memview.shape[0])

    return arr
