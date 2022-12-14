import numpy as np
import pandas as pd
import cmath

from FDCLF.FDCLF_functions import Ybus_fdclf, iterate_fdclf, printing_B_dash, printing_B_double_dash, Update_V_vec
from Newton_raphson.NR_functions import read_buses, Ybus, printing_buses, printing_Y_bus, Q_violated, Q_max_violation, printing_lines, PQ_to_PV  
from Newton_raphson.NR_network import Network


bus_vec = read_buses('PartA/Busdata.csv')
power_network = Network(bus_vec)
Ybus = Ybus('PartA/impedances.csv', len(bus_vec))


convergence = 0.00001
Q_max = [0.5, 5, 1.5,5,5]


def FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method):
    V_vec_1, V_vec_2 = power_network.get_V_vec_FD() #Returns V_vec_1 to calculate Q, and V_vec_2 to caluculate P. 
    Q_vec_FD = power_network.get_Q_vec_FD() #Returns Q_vec without slack Q
    P_vec_FD = power_network.get_P_vec_FD() #Returns P_vec without Slcak and PV P. 
    bus_num_init = power_network.get_bus_num_vec()  #Returns bus number as a vector
    num_buses = len(bus_num_init) #Calculates how many busses there are in the network. 
    delta_vec_init = power_network.get_delta_vec_FD() #Returns delta_values. without slack
    bus_type_vec = power_network.get_bus_type_vec() #Returns vector with bus type. O if slack, 1 if PV and 2 if PQ.
    bus_type_init_clean = power_network.get_bus_type_vec() #Returns vector with bus type. O if slack, 1 if PV and 2 if PQ. No change
    V = power_network.get_V_calc() #Returns V_values. If unknown, returns 1
    delta = power_network.get_delta_vec() #Returns delta_values. If unknown, returns 0
    delta = np.zeros(len(bus_vec)) #For Q_max after 
    B_dash, B_double_dash = Ybus_fdclf(bus_type_vec, len(bus_vec), Ybus) #Returns the FDLF Y-matrixes
    printing_Y_bus(Ybus)
    printing_B_dash(B_dash) 
    printing_B_double_dash(B_double_dash)

    delta_Delta = [1]
    delta_V = [1]
    i= 0
    while((abs(max(np.real(delta_Delta))) > convergence) or (abs(max(np.real(delta_V))) > convergence)):
        if (i==0):
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_Delta, delta_V, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network, bus_type_vec, Q_vec_FD, P_vec_FD = iterate_fdclf(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec_init, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, Q_max, power_network, method, Q_limit, reactive_limits_method)
                                                                             
            i += 1
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
        
        else:
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_Delta, delta_V, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network, bus_type_vec, Q_vec_FD, P_vec_FD = iterate_fdclf(num_buses, bus_num_init, V_updated, V_vec_1_updated, V_vec_2_updated, delta_updated, delta_updated, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, Q_max, power_network, method, Q_limit, reactive_limits_method)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            
            i += 1
    #Her testes Q_violation etter endt iterasjon. Da m?? ikke Q_violated v??re aktivert inne i while-loopen over. 
    if(reactive_limits_method == 'after'):    
        if(Q_violated(Q_max, Q_updated, bus_type_vec)):
                Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
                power_network = PQ_to_PV(bus_type_init_clean, bus_type_vec, power_network, V_updated)
                bus_type_vec = power_network.get_bus_type_vec()
                V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
                Q_vec_FD = power_network.get_Q_vec_FD()
                P_vec_FD = power_network.get_P_vec_FD()
                V_vec_1_updated, V_vec_2_updated = Update_V_vec(bus_type_vec, V_vec_1, V_vec_2, V_updated)
        delta_Delta = [1]
        delta_V = [1]
        while((abs(max(np.real(delta_Delta))) > convergence) or (abs(np.real(max(delta_V))) > convergence)):
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_Delta, delta_V, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network, bus_type_vec, Q_vec_FD, P_vec_FD = iterate_fdclf(num_buses, bus_num_init, V_updated, V_vec_1_updated, V_vec_2_updated, delta_updated, delta_updated, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, Q_max, power_network, method, Q_limit, reactive_limits_method)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            i += 1
    Power_network = PQ_to_PV(bus_type_init_clean, bus_type_vec, power_network, V_updated) #Sets the transfrormed PV_bus back to a PV_bus.
    printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_init_clean)
    printing_lines(bus_vec, "PartA/impedances.csv", V_updated, Ybus, delta_updated)
    
    return P_updated, Q_updated




