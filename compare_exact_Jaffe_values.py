#!/usr/bin/python
#

import sys
import os
import numpy as np
import math
import csv
import optparse
from math import sqrt
from Jaffe_Diode import Jaffe_Diode 

cmdp = optparse.OptionParser()
cmdp.add_option("--inputFile", action="store", type="string", dest="inputFile")
#cmdp.add_option("--outputFile", action="store", type="string", dest="outputFile")

opt, args = cmdp.parse_args()

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

# print( "Parameters: ", "d  = ", d, " Jfac = ", Jfac, " beam_energy (eV) = ", beam_energy, " Mach = ", Mach0, " V1 = ", V1)
# print( "Nx = ", Nx, "Ny = ", Ny, " Nz = ", Nz)
# print( "J = ", J, "n0 = ", n0, "u0 = ", u0)
# print( "phi0 = ", phi0, "E0 = ", E0)
# print( "P0 = ", P0, "T0 = ", T0, " cs0 = ", cs0, " tauf0 = ", tauf0)
# print( "J0 = ", J0)

E_err_L2_squared = 0.
E_L2_squared     = 0.
E_err_L1         = 0.
E_L1             = 0.
E_err_infty      = 0.
E_infty          = 0.
phi_err_L2_squared = 0.
phi_L2_squared     = 0.
phi_err_L1         = 0.
phi_L1             = 0.
phi_err_infty      = 0.
phi_infty          = 0.
ne_err_L2_squared = 0.
ne_L2_squared     = 0.
ne_err_L1         = 0.
ne_err_infty      = 0.
ne_infty          = 0.
ne_L1             = 0.
vx_err_L2_squared = 0.
vx_L2_squared     = 0.
vx_err_L1         = 0.
vx_L1             = 0.
vx_err_infty      = 0.
vx_infty          = 0.
domain_size      = 0.

# Write outputFile with x values from input
with open(opt.inputFile, 'r') as f:
#    with open(opt.outputFile, "w") as outputFile:
        reader = csv.reader(f, delimiter='\t')
        reader.next() # skip header
        for row in reader:
            x      = float(row[0])
            weight = float(row[1])
            E_actual   = float(row[2])
            phi_actual = float(row[3])
            ne_actual  = float(row[4])
            vx_actual  = float(row[5])
            phi_exact  = jaffe.compute_phi(x)
            ne_exact   = jaffe.compute_ne(x)
            E_exact    = jaffe.compute_E(x)
            vx_exact   = jaffe.compute_velocity(x)
            E_err_L2_squared   += (E_exact - E_actual) ** 2 * weight
            E_L2_squared       += E_exact **2 * weight
            phi_err_L2_squared += (phi_exact - phi_actual) ** 2 * weight
            phi_L2_squared     += phi_exact **2 * weight
            ne_err_L2_squared  += (ne_exact - ne_actual) ** 2 * weight
            ne_L2_squared      += ne_exact **2 * weight
            vx_err_L2_squared  += (vx_exact - vx_actual) ** 2 * weight
            vx_L2_squared      += vx_exact **2 * weight
            E_err_L1           += abs(E_exact - E_actual) * weight
            E_L1               += abs(E_exact) * weight
            phi_err_L1         += abs(phi_exact - phi_actual) * weight
            phi_L1             += abs(phi_exact) * weight
            ne_err_L1          += abs(ne_exact - ne_actual) * weight
            ne_L1              += abs(ne_exact) * weight
            vx_err_L1          += abs(vx_exact - vx_actual) * weight
            vx_L1              += abs(vx_exact) * weight
            E_err_infty   = max(E_err_infty  , abs(E_exact - E_actual))
            E_infty       = max(E_infty      , abs(E_exact))
            phi_err_infty = max(phi_err_infty, abs(phi_exact - phi_actual))
            phi_infty     = max(phi_infty    , abs(phi_exact))
            ne_err_infty  = max(ne_err_infty , abs(ne_exact - ne_actual))
            ne_infty      = max(ne_infty     , abs(ne_exact))
            vx_err_infty  = max(vx_err_infty , abs(vx_exact - vx_actual))
            vx_infty      = max(vx_infty     , abs(vx_exact))
            domain_size        += weight
 #           print >>outputFile, '\t'.join([str(x),str(E_exact),str(phi_exact),str(ne_exact),str(vx_exact)])

E_abs_err = sqrt(E_err_L2_squared)
E_rel_err = E_abs_err / sqrt(E_L2_squared)
phi_abs_err = sqrt(phi_err_L2_squared)
phi_rel_err = phi_abs_err / sqrt(phi_L2_squared)
ne_abs_err = sqrt(ne_err_L2_squared)
ne_rel_err = ne_abs_err / sqrt(ne_L2_squared)
vx_abs_err = sqrt(vx_err_L2_squared)
vx_rel_err = vx_abs_err / sqrt(vx_L2_squared)

print("domain size %10.3E" % domain_size)
print("")
print("L^2 norms:")
print("E:       %10.3E" % sqrt(E_L2_squared))
print("phi:     %10.3E" % sqrt(phi_L2_squared))
print("ne:      %10.3E" % sqrt(ne_L2_squared))
print("vx:      %10.3E" % sqrt(vx_L2_squared))
print("")
print("absolute L^2 error:")
print("E:       %10.3E" % E_abs_err)
print("phi:     %10.3E" % phi_abs_err)
print("ne:      %10.3E" % ne_abs_err)
print("vx:      %10.3E" % vx_abs_err)
print("")
print("relative L^2 error:")
print("E:       %10.3E" % E_rel_err)
print("phi:     %10.3E" % phi_rel_err)
print("ne:      %10.3E" % ne_rel_err)
print("vx:      %10.3E" % vx_rel_err)
print("")
print("")
print("absolute L^1 error:")
print("E:       %10.3E" % E_err_L1)
print("phi:     %10.3E" % phi_err_L1)
print("ne:      %10.3E" % ne_err_L1)
print("vx:      %10.3E" % vx_err_L1)
print("")
print("relative L^1 error:")
print("E:       %10.3E" % (E_err_L1 / E_L1))
print("phi:     %10.3E" % (phi_err_L1 / phi_L1))
print("ne:      %10.3E" % (ne_err_L1 / ne_L1))
print("vx:      %10.3E" % (vx_err_L1 / vx_L1))
print("")
print("")
print("absolute L^infty error:")
print("E:       %10.3E" % E_err_infty)
print("phi:     %10.3E" % phi_err_infty)
print("ne:      %10.3E" % ne_err_infty)
print("vx:      %10.3E" % vx_err_infty)
print("")
print("relative L^infty error:")
print("E:       %10.3E" % (E_err_infty / E_infty))
print("phi:     %10.3E" % (phi_err_infty / phi_infty))
print("ne:      %10.3E" % (ne_err_infty / ne_infty))
print("vx:      %10.3E" % (vx_err_infty / vx_infty))

# Compute minimum potential

# phimin_exact = jaffe.compute_phimin()
# print("phimin_exact = ", phimin_exact)
# print(" ")
# 
# # Estimate particle transit time
# tau = 0.0
# #cfl = 0.0
# #cflmin =  1.0E10
# #cflmax = -1.0E10
# vxmin  =  1.0E10
# vxmax  = -1.0E10
# 
# Nnd = Nnode    #  For accuracy of cflmin, cflmax 
# dxt = d/(Nnd-1)
# 
# for i in range(Nnd-1):
#      xa  = i*dxt
#      xb  = (i+1)*dxt
#      ua  = jaffe.compute_velocity(xa)
#      ub  = jaffe.compute_velocity(xb)
# #    tau = tau + (xb-xa)/( 0.5*(ua+ub) )
#      tau = tau + (xb-xa)*0.5*(ua+ub)/(ua*ub)
#      # cfl must use mesh resolution not Jaffe integration resolution
# #    cfl = dt*(0.5*(ua+ub))/(xb-xa)
# #    if (cfl < cflmin): cflmin = cfl
# #    if (cfl > cflmax): cflmax = cfl
#      if (ua  < vxmin): vxmin = ua
#      if (ua  > vxmax): vxmax = ua
# 
# #print("tau,cflmin,cflmax,vxmin,vxmax: ",tau, " ", cflmin, " ", cflmax, " ", vxmin, " ", vxmax)
# print("tau,vxmin,vxmax: ",tau, " ", vxmin, " ", vxmax)
# 
# # Fini!
