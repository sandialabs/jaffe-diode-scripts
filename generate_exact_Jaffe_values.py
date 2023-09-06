# 
#!/usr/bin/env python
#

import sys
import os
import numpy as np
import math
from Jaffe_Diode import Jaffe_Diode 

def getopts(argv):
    opts = {}  # Empty dictionary to store key-value pairs.
    while argv:  # While there are arguments left to parse...
        argc = len(argv)
        if argv[0][0] == '-' and argc > 1:  # Found a "-name value" pair.
            opts[argv[0]] = argv[1]  # Add key and value to the dictionary.
        else:
            opts[argv[0]] = ' '  # Add key and value to the dictionary.
        argv = argv[1:]  # Reduce the argument list by copying it starting from index 1.
    return opts

from sys import argv

################
# Parsing 
################
#print(sys.argv) 
#argc = len(sys.argv)
#print ('argc =', argc)

#if  argv.count('-h') >= 1:
#    print('\nusage: python diode_driver.py -t time -Nt num_timesteps')
#    exit(1)

#myargs = getopts(argv)
##print myargs

#if '-t' not in myargs:  
#    print('\nusage: python diode_driver -t time -Nt num_timesteps')
#    exit(1)

#if '-Nt' not in myargs:  
#    print('\nusage: python diode_driver.py -t time -Nt num_timesteps')
#    exit(1)

#time = float(myargs['-t']) # (sec)
#Nt   = int(myargs['-Nt'])  
#dt   = time/Nt 
#print "time = ", time, " Nt = ", Nt, " dt = ", dt

# Some Units
# 1 Joule                = 1 kg.m^2/sec^2
# 1 Volt                 = 1 kg.m^2/(Amp.s^3) = 1 Joule/Coulomb
# 1 Coulomb              = 1 Ampere.sec (Amp.s)
# 1 Ampere               = 1 Coulomb/s (C/s)
# 1 Volt                 = 1 kg.m^2/(Amp.s^3) = 1 Joule/Coulomb
# electrical charge (qe) = 1 Coulomb = 1 Amp.s
# electron volt (eV)     = 1.6021766208E-19 Joules

# EMPIRE-Fluid Parameters:
qe   = 1.60217662e-19  # elementary charge = (e = 1.60217662e-19) C
me   = 9.10938356e-31  # electron mass     = (m_e = 9.10938356e-31) kG
eps0 = 8.85418781e-12  # permeability of free space in SI
kB   = 1.38064852e-23  # Boltzmann's constant
adiabatic_index = 1.01 # gamma

Nnode = 101
Nx   = 100
Ny   = 0
Nz   = 0
d    = 0.01            # (m)
dx   = d/(Nnode-1)     # (m)
xc0  = dx/2.0
Ly   = 0.0             # (m)
Lz   = 0.0             # (m)

# Case 1
#V1         = 100.0      # (Volts)
#Jfac        = 1.0       # nond
#beam_energy = 1.0*qe*V1 # (eVolts)
#Mach0       = 1.25

# Case 2
Jfac = 1.0             # nond
beam_energy = 1.0E4*qe # (eVolts)
V1     = 0.0      # (Volts)
Mach0  = 1.25

J0    = (16.0/9.0)*eps0*math.sqrt( 2.0*qe/me )*math.pow(beam_energy/qe,1.5)/(d*d)
J     = Jfac*J0    # (Amps/m^2/s)

jaffe = Jaffe_Diode(d, J, beam_energy, V1)

u0    = math.sqrt(2.0*beam_energy*qe/me)
n0    = J/(qe*u0)
E0    = jaffe.compute_E(xc0)
phi0  = jaffe.compute_phi(xc0)
T0    = u0*u0*me/(Mach0*Mach0*adiabatic_index*kB)
P0    = n0*kB*T0
cs0   = math.sqrt(adiabatic_index*kB*T0/me)
tauf0 = d/abs(u0)

print( "Parameters: ", "d  = ", d, " Jfac = ", Jfac, " beam_energy (eV) = ", beam_energy, " Mach = ", Mach0, " V1 = ", V1)
print( "Nx = ", Nx, "Ny = ", Ny, " Nz = ", Nz)
print( "J = ", J, "n0 = ", n0, "u0 = ", u0)
print( "phi0 = ", phi0, "E0 = ", E0)
print( "P0 = ", P0, "T0 = ", T0, " cs0 = ", cs0, " tauf0 = ", tauf0)
print( "J0 = ", J0)

# Write Jaffe.out Profiles
with open("jaffe.out", "w") as F1:
  for i in range(Nnode-1):
     xc  = xc0 + dx*i
     phi = jaffe.compute_phi(xc)
     ne  = jaffe.compute_ne(xc)
     E   = jaffe.compute_E(xc)
     vx  = jaffe.compute_velocity(xc)
#    print("i,xc,phi,E,ne,vx: ", i, " ", xc, " ", phi, " ", E, " ", ne, " ", vx)
     
     F1.write("{0:25.14e}{1:25.14e}{2:25.14e}{3:25.14e}{4:25.14e}\n"
       .format(xc, phi, E, ne, vx))

# Compute minimum potential

phimin_exact = jaffe.compute_phimin()
print("phimin_exact = ", phimin_exact)
print(" ")

# Estimate particle transit time
tau = 0.0
#cfl = 0.0
#cflmin =  1.0E10
#cflmax = -1.0E10
vxmin  =  1.0E10
vxmax  = -1.0E10

Nnd = Nnode    #  For accuracy of cflmin, cflmax 
dxt = d/(Nnd-1)

for i in range(Nnd-1):
     xa  = i*dxt
     xb  = (i+1)*dxt
     ua  = jaffe.compute_velocity(xa)
     ub  = jaffe.compute_velocity(xb)
#    tau = tau + (xb-xa)/( 0.5*(ua+ub) )
     tau = tau + (xb-xa)*0.5*(ua+ub)/(ua*ub)
     # cfl must use mesh resolution not Jaffe integration resolution
#    cfl = dt*(0.5*(ua+ub))/(xb-xa)
#    if (cfl < cflmin): cflmin = cfl
#    if (cfl > cflmax): cflmax = cfl
     if (ua  < vxmin): vxmin = ua
     if (ua  > vxmax): vxmax = ua

#print("tau,cflmin,cflmax,vxmin,vxmax: ",tau, " ", cflmin, " ", cflmax, " ", vxmin, " ", vxmax)
print("tau,vxmin,vxmax: ",tau, " ", vxmin, " ", vxmax)

# Fini!
