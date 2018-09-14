import numpy as np
cimport numpy as np

## References the functions defined in the header of the C library
cdef extern from "multiply.h":
    void multiply_by_10_in_C(double arr[], unsigned int n)


def py_multiply_by_10_in_C(arr) -> None:
    if not arr.flags['C_CONTIGUOUS']:
        arr = np.ascontiguousarray(arr) 

    ## Creates a contiguous view (ie. pointer) to the contiguous array
    cdef double[::1] arr_memview = arr

    multiply_by_10_in_C(&arr_memview[0], arr_memview.shape[0])

    return arr
