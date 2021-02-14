"""
Sample script to compute ignition delay times of fuel-air mixtures for
constant-pressure and constant-volume homogeneous reactors

Tested with python2 and Cantera 2.3

Author: Francisco E. Hernandez Perez
"""

###################################################################
#  Load modules
###################################################################
import cantera as ct
import numpy as np
import matplotlib.pyplot as plt


###################################################################
# Define some functions
###################################################################

##################################################
# Generate list of values
##################################################
def Generate_list(min_val, max_val, delta):
    n_pts = (max_val-min_val)/delta + 1
    tmp_list = []
    for i in range( int(n_pts) ):
        tmp_list.append( min_val + i*delta )
    if tmp_list[-1] < max_val:
        tmp_list.append( max_val )    
    return tmp_list


##################################################
# Mixture composition from equivalence ratio
##################################################
def composition_CxHyOz_Air(gas, fuel_species, phi):
    iFuel = gas.species_index(fuel_species)
    iO2 = gas.species_index('O2')
    iN2 = gas.species_index('N2')
    C_atoms = 0.0
    H_atoms = 0.0
    O_atoms = 0.0
    if 'C' in gas.element_names: C_atoms = gas.n_atoms(fuel_species, 'C')
    if 'H' in gas.element_names: H_atoms = gas.n_atoms(fuel_species, 'H')
    if 'O' in gas.element_names: O_atoms = gas.n_atoms(fuel_species, 'O')
    stoich_O2 = C_atoms + H_atoms/4.0 - O_atoms/2.0
    #print('stoich_O2 = %g'%stoich_O2)
    X = np.zeros(gas.n_species)
    X[iFuel] = phi
    X[iO2] = stoich_O2
    X[iN2] = 3.76*stoich_O2
    #print 'Composition:',fuel_species,'=',X[iFuel],', O2=',X[iO2],', N2=',X[iN2]
    return X


##################################################
# Heat release rate [W/m^3]
##################################################
def heat_release_rate(reactor):
    wdot_molar = reactor.kinetics.net_production_rates # [kmol/m^3/s]
    h_molar = reactor.thermo.partial_molar_enthalpies  # [J/kmol]
    nsp = reactor.thermo.n_species
    hrr = 0.0
    for i in range(nsp):
        hrr += -(wdot_molar[i]*h_molar[i]) 
    ###########
    return hrr


##################################################
# Reactor time-stepping
##################################################
def RunReactor(gas_filename, const_pressure, fuel_species, temp, press, phi, time_step, time_max):
    # set the gas with its initial condictions and composition
    gas = ct.Solution(gas_filename)
    gas.TPX = temp, press, composition_CxHyOz_Air(gas, fuel_species, phi)


    # determine equilibrium temperature, which will be used later to stop the simulation
    gas2 = ct.Solution(gas_filename)
    gas2.TPX = temp, press, composition_CxHyOz_Air(gas, fuel_species, phi)
    if const_pressure == True:
        gas2.equilibrate('HP')
    else:
        gas2.equilibrate('UV')

    
    # set the reactor
    if const_pressure == True:
        r = ct.IdealGasConstPressureReactor(gas)  # constant pressure
    else:
        r = ct.IdealGasReactor(gas)  # constant volume            

    # set the simulator    
    sim = ct.ReactorNet([r])

    # initialize some parameters
    Teq = gas2.T
    Told = r.T
    dT = 0.0
    dTmax = 1.0
    HHRmax = 0.0

    time = 0.0
    timeig_dT400 = 0.0
    time_maxdT = 0.0
    time_max_HRR = 0.0
    ignited = False

    # determine output file name
    stem, ext = gas_filename.split('.')
    out_file = stem + '-P_' + str(press/ct.one_atm) + '-phi_' + str(phi) + '-T_' + str(temp) + '.dat'
    
    file = open(out_file,'w')
    # header info for the output file
    file.write('#t(s) T(K) P(atm) HRR(W/m3)')
    for j in range(gas.n_species): file.write(' %s'%(gas.species_names[j]))
    file.write('\n')

    # time-march the system of ODEs and determine the ignition delay
    # time based on three criteria: T=To+400K (timeig_dT400), maximum rate of change of
    # temperature (time_maxdT) and peak of heat release rate (time_max_HRR)
    while True:
        time += time_step
        sim.advance(time)
        # compute heat release rate
        hrr = heat_release_rate(r)
        # check for max value of heat release
        if hrr > HHRmax:
            HHRmax = hrr
            time_max_HRR = time
        dT = r.T-Told
        if dT > dTmax:
            dTmax = dT
            time_maxdT = time
        Told = r.T
        if ignited == False and r.T > (temp+400.0):
            timeig_dT400 = time
            ignited = True
        # output to file
        file.write('%10.3e %10.3f %10.3e %10.3e'%(sim.time, r.T, r.thermo.P/ct.one_atm, hrr))
        for j in range(gas.n_species): file.write(' %10.3e'%(gas.Y[j]))
        file.write('\n')
        if (time > time_max) or (r.T >= (Teq-5.0)):
            break 
    file.close
    #print('Ignition time: %g'%timeig)
    return timeig_dT400, time_maxdT, time_max_HRR 


###################################################################
# Main
###################################################################
if __name__ == "__main__":
    # define the gas, fuel and reactor type
    gas_filename = "h2_burke.xml"
    fuel_species =  'H2'
    const_pressure = True  # True or False
    
    # set time-stepping parameters
    time_max =  10.0
    time_step = 1e-5

    # set range of pressure values in atm
    min_val = 10.0
    max_val = 10.0 
    step = 1.0
    P_values = Generate_list(min_val, max_val, step)
    # convert to Pascals
    for i in range(len(P_values)):
        P_values[i] *= ct.one_atm
    print '\nInitial pressure values:',P_values
    
    # set range of equivalence ratio values
    min_val = 1.0 
    max_val = 1.0 
    step = 0.05
    phi_values = Generate_list(min_val, max_val, step)
    print '\nEquivalence ratio values:',phi_values

    # set range of temperature values
    min_val = 800.0
    max_val = 1200.0 
    step = 20.0
    T_values = Generate_list(min_val, max_val, step)
    print '\nInitial temperature values:',T_values
    print '\n'

    n_Tvals = len(T_values)
    ignition_time = np.zeros(n_Tvals)
    # main loop
    fig1 = plt.figure(1)
    # loop over pressure
    for P in P_values:
        # loop over equivalence ratio
        for phi in phi_values:
            # loop over temperature
            for i in range(n_Tvals):
                # run the reactor
                timeig_dT400, time_maxdT, time_max_HRR = RunReactor(gas_filename, const_pressure, fuel_species, T_values[i], P, phi, time_step, time_max)
                print 'Ignition time for P =',P, ', phi =',phi, ', T =',T_values[i], '==>', timeig_dT400, ',', time_maxdT, ',', time_max_HRR, 's'
                time_ig = timeig_dT400
                ignition_time[i] = time_ig
            print '\nignition_delay:',ignition_time
            # plots
            plt.figure(1)
            p1, = plt.plot(T_values,ignition_time)
            plt.xlabel('Temperature [K]',fontsize=14)
            plt.ylabel('t_ig [s]',fontsize=14)
            plt.grid(True)
            str_leg = 'phi=' + str(phi)
            #plt.legend([p1], [str_leg], loc=1)
            
    # show plots
    plt.show()

