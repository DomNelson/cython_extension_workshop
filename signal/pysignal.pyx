import numpy as np
cimport numpy as np


def sort_ped(Ped):
    ninds = len(Ped.inds)
    sorted_ped_arr = np.zeros((ninds, 3), dtype=np.int32)

    for i in range(ninds):
        ind = Ped.inds[i]
        father = Ped.fathers[i]
        mother = Ped.mothers[i]

        father_ix = -1
        mother_ix = -1
        if father != 0:
            father_ix = Ped.ind_dict[father]
        if mother != 0:
            mother_ix = Ped.ind_dict[mother]

        sorted_ped_arr[i] = [ind, father_ix, mother_ix]

    return sorted_ped_arr


## References the functions defined in the header of the C library
cdef extern from "signal.h":
    ## This sets the attributes accessible in the class attribute
    ## of type node_t below
    ctypedef struct node_t:
        int ID
        float weight

    ctypedef struct ped_t:
        pass

    ## TODO: Double check int/uint casting here
    void multiply_by_10_in_C(double arr[], unsigned int n)
    void print_node(node_t node) 
    ped_t *ped_alloc()
    int ped_nodes_alloc(ped_t *ped, int num_nodes, int num_samples)
    int free_ped(ped_t *ped)
    int ped_load(ped_t *ped, int *inds, int *fathers, int *mothers, int num_inds)
    int ped_print_nodes(ped_t *ped)
    int ped_samples_alloc(ped_t *ped, int num_samples)
    int ped_load_samples_from_idx(ped_t *ped, int *samples_idx, int *genotypes, int num_samples)
    int ped_init_sample_weights(ped_t *ped)
    int ped_print_samples(ped_t *ped)

    int ped_set_all_weights(ped_t *ped, double val)
    double ped_get_node_weight_from_idx(ped_t *ped, int node_idx)
    int update_parent_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx)
    int update_parent_not_carrier_from_idx(ped_t *ped, int node_idx, int sample_idx)

    int ped_climb_step(ped_t *ped)


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


## Classes best used when we need a persistent object to
## pass to functions/manipulate, and which uses a custom
## type that can't simply be referenced in Python
cdef class cPed:
    cdef ped_t *ped
    cdef int ret

    def __cinit__(self):
        self.ped = ped_alloc()


    ## For this to work, the memory for the ped_t self.ped struct
    ## had to be allocated with calloc. Then it could be freed, after which
    ## self.ped would properly evaluate as NULL
    def __dealloc__(self):
        if self.ped is not NULL:
            free_ped(self.ped)


    def load_samples(self, sample_arr, genotypes_arr=None):
        num_samples = len(sample_arr)
        ret = ped_samples_alloc(self.ped, num_samples)
        if ret != 0:
            raise MemoryError()

        if genotypes_arr is None:
            print "Setting all samples as heterozygotes"
            genotypes_arr = np.ones(num_samples, dtype=np.int32)

        cdef int [::1] samples = np.ascontiguousarray(sample_arr)
        cdef int [::1] genotypes = np.ascontiguousarray(genotypes_arr)
        ret = ped_load_samples_from_idx(self.ped, &samples[0],
                &genotypes[0], num_samples)
        if ret != 0:
            ## Using this as generic error
            raise MemoryError()


    def load_ped(self, ped_arr, num_samples):
        num_nodes = len(ped_arr)
        ret = ped_nodes_alloc(self.ped, num_nodes, num_samples)
        if ret != 0:
            raise MemoryError()

        ## Should implement a check that the array is properly sorted
        if not ped_arr.flags['C_CONTIGUOUS']:
            ## Makes a contiguous copy of the numpy array
            ped_arr = np.ascontiguousarray(ped_arr) 

        ## - Should better explain what's happening here - takes a numpy 'view'
        ## of the array I think, using special syntax to ensure it is
        ## contiguous
        ##  - Note also that passing arrays of pointers is not straightforward,
        ## (numpy ndarrays are always 'flat') so we split into multiple 1D
        ## arrays here
        cdef int [::1] inds = np.ascontiguousarray(ped_arr[:, 0])
        cdef int [::1] fathers = np.ascontiguousarray(ped_arr[:, 1])
        cdef int [::1] mothers = np.ascontiguousarray(ped_arr[:, 2])

        ret = ped_load(self.ped, &inds[0], &fathers[0], &mothers[0], ped_arr.shape[0])
        if ret != 0:
            ## Using this as generic error
            raise MemoryError()

        if num_nodes < 50:
            for i, row in enumerate(ped_arr):
                print i, row
        else:
            print num_nodes, "individuals loaded"

    ## Defining with cpdef (also regular def, which has more overhead) exposes
    ## the function to the Python API
    cpdef print_nodes(self):
        ped_print_nodes(self.ped)

    cpdef print_samples(self):
        ped_print_samples(self.ped)

    cpdef set_all_weights(self, val):
        ped_set_all_weights(self.ped, val)

    cpdef get_node_weight(self, ind_idx):
        return ped_get_node_weight_from_idx(self.ped, ind_idx)

    cpdef update_parent_carrier(self, ind_idx, sample_idx):
        update_parent_carrier_from_idx(self.ped, ind_idx, sample_idx)

    cpdef update_parent_not_carrier(self, ind_idx, sample_idx):
        update_parent_not_carrier_from_idx(self.ped, ind_idx, sample_idx)

    cpdef climb_step(self):
        ped_climb_step(self.ped)

    cpdef init_sample_weights(self):
        ped_init_sample_weights(self.ped)
