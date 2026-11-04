import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import sys

N = 18000
if len(sys.argv) == 2:
    N = int(sys.argv[1])
D = 1
totalsteps = 200000
if len(sys.argv) == 3:
    totalsteps = int(sys.argv[2])

cols = np.arange(6, 6+N**D, 1)
def get_array(filename, **kwargs):
    return np.genfromtxt( filename, usecols=cols, dtype=complex, **kwargs )

count = 6 + N
def get_melement(content):
    arr = content.split(" ")
    data = np.array(arr[6:]).astype(complex)
    return data

def get_matrix(filename, groundstate, header=23):
    matrix_element = np.empty(totalsteps, dtype=complex)
    rows = np.arange(23, 23 + totalsteps, 1)
    x = np.arange(-N/2, N/2, 1)

    gx = np.conj(groundstate) * x
    
    with open(filename) as file:
        # skip header
        for i in range(23):
            next(file)
        
        # get matrix elements
        for i, line in enumerate(file):
            # stop iteration early if an error occurs while reading data
            try:
                data = get_melement(line)
                matrix_element[i] = np.dot(gx, data)
            except Exception as err:
                import traceback
                print(f"Error occured, iteration stopped early at row {i}")
                traceback.print_exc()
                break # or raise StopIteration?

    matrix_element[i:] = np.nan
    return matrix_element

gstate = get_array(f"data/gstate_N{N}_epssq1e-7.txt", max_rows=1)
melement = get_matrix(f"data/N{N}_epssq1e-7.txt", gstate)
np.savetxt(f"melementN{N}.txt", melement, header=f"N={N}\neps2=1e-7\ntau=5e-5\ntotal=200000")

