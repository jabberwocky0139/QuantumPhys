# coding: utf-8

import os
import re
import numpy as np
from scipy.integrate import quad
from scipy.fftpack import fft
from scipy.fftpack import ifft
from inspect import getsource


# --*-- Set functions --*--

# Set potential
# Note1 : Active line must be attached semicolon(;) at the end of line
# Note2 : Use the math function which work on gnuplot as well without any change
def Potential(x):
    return 0.01*(x**2 - 8*abs(x) + 16);

# Potential on exp
def V(x):
    return -0.5j*Potential(x)

# Initial function
def Psi_0(x):
    y = lambda z: np.exp(-2*(z-7)**2)
    return y(x)/quad(y, -np.inf, np.inf)[0]/2

# file writer
def file_writer(filename, arr_func):
    with open(filename, "w") as f:
        for n, x in enumerate(np.linspace(-L/2, L/2, N)):
            print("{0}\t{1}".format(x, abs(arr_func[n])), file=f)


# --*-- Set constants and variables for initial condition --*--

# Space step, Volume, Time step
N, L, dt = 512, 25.0, 0.5

# Set expK
expK = np.array([])
for n in range(0, N):
    if(n < N/2):
        expK = np.append(expK, np.exp(-0.5j*(2*np.pi*n/N)**2*dt))
    else:
        expK = np.append(expK, np.exp(-0.5j*(2*np.pi*(n-N)/N)**2*dt))

# Set expV, initial function arr_Psi_old
expV ,arr_Psi = np.array([]), np.array([], dtype=complex)
for x in np.linspace(-L/2, L/2, N):
    expV = np.append(expV, np.exp(V(x)*dt))
    arr_Psi = np.append(arr_Psi, Psi_0(x))


# --*-- Time propagation by symplectic numerical solution --*--

# Time length, Output interbal, Total file piece
times, total_file = 40000, 500
M = times//total_file
for i in range(times):
    # output time depend function
    if(i%M == 0):
        file_writer("output{0}.txt".format(i//M), arr_Psi)
    
    # Multipling by time propagator of potential term
    arr_Psi *= expV
    
    # Fourier transformation
    arr_Psi = fft(arr_Psi)
    
    # Multipling time propagator of kinetic term
    arr_Psi *= expK
    
    # Inverse fourier transformation
    arr_Psi = ifft(arr_Psi)


# --*-- Gnuplot --*--

# Get potential function as string
pattern = re.compile("(return.+;)")
m = pattern.search(getsource(Potential))
str_potential = str(m.group(0))[7:-1]

# System call for gnuplot
gnuplot_call = 'gnuplot -e '
start = '"'
set_range = 'set xr[{0}:{1}]; set yr[{2}:{3}]; '.format(-L/2, L/2, 0, 0.5)
plot_initial = 'plot \\"output0.txt\\" w l lw 2 title \\"t = 0\\",{0} title \\"potential\\" w l lw 2; pause 2; '.format(str_potential)
do_for_declaration = 'do for[i = {0}:{1}:1]'.format(1, total_file-1)
do_for_start = "{"
do_for_procedure = 'plot sprintf(\\"output%d.txt\\", i) title sprintf(\\"t = %d\\", i) w l lw 2, {0} title \\"potential\\" w l lw 2; pause 0.05;'.format(str_potential)
end = '}"'

os.system("".join([gnuplot_call, start, set_range, plot_initial, do_for_declaration, do_for_start, do_for_procedure, end]))

os.system("rm output*.txt")

