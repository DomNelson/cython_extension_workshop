import numpy as np
import pymultiply


def multiply_by_10_in_python(arr):
    for i in range(len(arr)):
        arr[i] *= 10


for_python = np.arange(10)
for_C = np.arange(10, dtype=float)

print(for_python)
print(for_C)

multiply_by_10_in_python(for_python)
pymultiply.py_multiply_by_10_in_C(for_C)

print(for_python)
print(for_C)
