import numpy as np
import pandas as pd
import cmath

from FDCLF.FDCLF_functions import Ybus_fdclf, iterate_fdclf, printing_B_dash, printing_B_double_dash, Update_V_vec
from Task4.Task4_NR_func import Ybus_trans
from Task4.Trans_network import Network
#from Newton_raphson.NR_functions import read_buses, printing_buses, printing_Y_bus, printing_lines
from Newton_raphson.NR_functions import read_buses, printing_buses, printing_Y_bus, printing_lines, Q_violated, Q_max_violation


bus_vec = read_buses('PartA/Task4/BusdataWith7Buses.csv')
power_network_trans = Network(bus_vec)
Ybus = Ybus_trans('PartA/Task4/impedancesPart4.csv', len(bus_vec))


convergence = 0.00001
Q_max = [0.5, 5, 1.5,5,5]


def FDCLF_trans(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method):
    V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
    Q_vec_FD = power_network.get_Q_vec_FD()
    P_vec_FD = power_network.get_P_vec_FD()
    bus_num_init = power_network.get_bus_num_vec()
    num_buses = len(bus_num_init)
    delta_vec_init = power_network.get_delta_vec_FD()
    bus_type_vec = power_network.get_bus_type_vec()
    V = power_network.get_V_calc()
    delta = power_network.get_delta_vec()
    delta = np.zeros(len(bus_vec))
    printing_Y_bus(Ybus) 

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
    #Her testes Q_violation etter endt iterasjon. Da må ikke Q_violated være aktivert inne i while-loopen over. 
    if(reactive_limits_method == 'after'):    
        if(Q_violated(Q_max, Q_updated, bus_type_vec)):
                Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
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
    
    printing_lines(bus_vec, 'PartA/Task4/impedancesPart4.csv', V_updated, Ybus, delta_updated)
    
    return P_updated, Q_updated


"""
def FDCLF_trans(Ybus, power_network, convergence, Q_max, Q_limit):
    V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
    Q_vec_FD = power_network.get_Q_vec_FD()
    P_vec_FD = power_network.get_P_vec_FD()
    b_dash, b_double_dash = Ybus_fdclf(bus_vec, len(bus_vec), Ybus)#Ybus_fdclf('PartA/Task4/impedancesPart4.csv', len(bus_vec), 'PartA/Task4/BusdataWith7Buses.csv', power_network, Ybus)
    bus_num_init = power_network.get_bus_num_vec()
    num_buses = len(bus_num_init)
    delta_vec_init = power_network.get_delta_vec_FD()
    bus_type_vec = power_network.get_bus_type_vec()
    V = power_network.get_V_calc()
    delta = power_network.get_delta_vec()
    delta = np.zeros(len(bus_vec))
    P = power_network.get_P_vec()
    Q = power_network.get_Q_vec()
    
    printing_Y_bus(Ybus)
    printing_B_dash(b_dash)
    printing_B_double_dash(b_double_dash)

    delta_Delta = [1]
    delta_P = [1]
    i= 0
    while((abs(max(np.real(delta_Delta))) > convergence) and (abs(max(np.real(delta_P))) > convergence)):
        if (i==0):
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_P, delta_Delta, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network = iterate_fdclf(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec_init, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q, Q_max, power_network)
                                                                             
            i += 1
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
        
        else:
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_P, delta_Delta, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network = iterate_fdclf(num_buses, bus_num_init, V_updated, V_vec_1_updated, V_vec_2_updated, delta_updated, delta_updated, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q, Q_max, power_network)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            
            i += 1
        
    return P_updated, Q_updated
"""


