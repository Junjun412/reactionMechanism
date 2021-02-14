"""
An opposed-flow ethane/air diffusion flame
"""

import csv
import cantera as ct
import numpy as np
#import matplotlib.pyplot as plt

# Input parameters
p = 1.0*ct.one_atm  # pressure
tin_f = 298.0  # fuel inlet temperature
tin_o = 310.0  # oxidizer inlet temperature
mdot_o = 0.145919  # kg/m^2/s  # 0.11963 m/s
mdot_f = 0.132676  # kg/m^2/s  # 0.12031 m/s

comp_o = 'O2:0.25, N2:0.75'  # air composition
comp_f = 'C2H4:0.9, N2:0.1'  # fuel composition

width = 0.01375 # Distance between inlets is 2 cm

loglevel = 1  # amount of diagnostic output (0 to 5)

# Create the gas object used to evaluate all thermodynamic, kinetic, and
# transport properties.
#gas = ct.Solution('sk99.xml', 'gri30_mix')
gas = ct.Solution('grimech30.xml')
gas.TP = gas.T, p

# Create an object representing the counterflow flame configuration,
# which consists of a fuel inlet on the left, the flow in the middle,
# and the oxidizer inlet on the right.
f = ct.CounterflowDiffusionFlame(gas, width=width)

# Set the state of the two inlets
f.fuel_inlet.mdot = mdot_f
f.fuel_inlet.X = comp_f
f.fuel_inlet.T = tin_f

f.oxidizer_inlet.mdot = mdot_o
f.oxidizer_inlet.X = comp_o
f.oxidizer_inlet.T = tin_o

# Set the boundary emissivities
f.set_boundary_emissivities(0.0, 0.0)
# Turn radiation off
f.radiation_enabled = False

f.set_refine_criteria(ratio=4, slope=0.05, curve=0.25, prune=0.04)

# Set transport model to Le=1 instead of mixture-averaged
#f.transport_model = 'UnityLewis'  #'Mix'

# Solve the problem
f.solve(loglevel, auto=True)
f.show_solution()
#f.mixture_fraction('C')
f.save('c2h4_diffusion.xml')
f.write_csv('c2h4_diffusion_norad.csv', species='Y', quiet=False)

#write the solution to a CSV file
#x = f.grid()
#Z = f.mixture_fraction('C')
T = f.T
#u = f.u
#V = f.V
#A1 = gas.massFraction('A1')
#fcsv = open('c2h4_diffusion.csv','w')
#writer = csv.writer(fcsv)
#writer.writerow(['Z (-)','u (m/s)','V (1/s)','T (K)'])
#for n in range(339):
#	f.setGasState(n)
#	writer.writerrow([Z[n], u[n], V[n], T[n] ]+list(gas.Y))
#	writer.writerow([Z[n], u[n], V[n], T[n] ])
#fcsv.close()
print ('Max Tem is {0:.2f} '.format(max(T)))
print ('solution saved to csv file')


# Define the element to follow in the reaction path diagram:
element = 'C'

# Initiate the reaction path diagram:
diagram = ct.ReactionPathDiagram(gas, element)

# Options for cantera:
diagram.show_details = False
diagram.font='CMU Serif Roman'
#diagram.threshold=0.01
diagram.dot_options='node[fontsize=20,shape="box"]'
diagram.title = 'Reaction path diagram following {0}'.format(element)

# Define the filenames:
dot_file = 'ReactionPathDiagram.dot'
#img_file = 'ReactionPathDiagram.png'

# Write the dot-file first, then create the image from the dot-file with customizable
# parameters:
diagram.write_dot(dot_file)


# write the velocity, temperature, and mole fractions to a CSV file
#f.write_csv('c2h6_diffusion.csv', quiet=False)

#f.show_stats(0)

# Plot Temperature without radiation
#figTemperatureModifiedFlame = plt.figure()
#plt.plot(f.flame.grid, f.T, label='Temperature without radiation')
#plt.title('Temperature of the flame')
#plt.ylim(0,2500)
#plt.xlim(0.000, 0.020)

# Turn on radiation and solve again
#f.radiation_enabled = True
#f.solve(loglevel=1, refine_grid=False)
#f.show_solution()

# Plot Temperature with radiation
#plt.plot(f.flame.grid, f.T, label='Temperature with radiation')
#plt.legend()
#plt.legend(loc=2)
#plt.savefig('./c2h6_diffusion.pdf')
