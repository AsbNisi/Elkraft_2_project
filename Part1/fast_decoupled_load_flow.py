import numpy as np
import pandas as pd
import cmath

from FDCLF_functions import Ybus_fdclf 
from NR_functions import read_buses, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV
from DCLF_functions import Ybus_dclf, printing_Y_bus, iterate_dclf
from NR_network import Network, Buses, PQ, VD


bus_vec = read_buses('Part1/Busdata.csv')
power_network = Network(bus_vec)
#print(power_network.get_bus(2))
Ybus = Ybus('Part1/impedances.csv', 5)

convergence = 0.00001
Q_max = [0.5, 5, 1.5,5,5]


def FDCLF(Ybus, power_network, convergence, Q_max):
    V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
    Q_vec_FD = power_network.get_Q_vec_FD()
    P_vec_FD = power_network.get_P_vec_FD()
    b_dash, b_double_dash = Ybus_fdclf('Part1/impedances.csv', 5, 'Part1/Busdata.csv', power_network)
    bus_num_init = power_network.get_bus_num_vec()
    num_buses = len(bus_num_init)
    delta_vec_init = power_network.get_delta_vec_FD()
    bus_type_vec = power_network.get_bus_type_vec()
    V = power_network.get_V_calc()

    
    delta_vd = [1,1,1,1,1,1,1]
    PQ_vec_updated = PQ_vec.copy()
    i= 0
    #print("V")
    #print(V)
    while(abs(max(np.real(delta_vd))) > convergence):
        if (i==0):
            print("Iteration", i+1, ": \n")
            bus_type = power_network.get_bus_type_vec()
            V_updated, delta_updated, P_updated, Q_updated = iterate_fdclf(num_buses, V, V_vec_1, V_vec_2, delta, delta_vec, V_vec, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash)
            i += 1
            #print("V_updated")
            #print(V_updated)
            #print("delta_updated")
            #print(delta_updated)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            """
            Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
            power_network = PQ_to_PV(bus_type_init, Q_updated, Q_max, V_updated, power_network)
            V = power_network.get_V_vec()
            Q_calc = power_network.get_Q_vec()
            VD_vec, VD_jacobian = power_network.get_VD_jacobian()
            PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
            delta = power_network.get_delta_calc()
            PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
            VD_vec, VD_jacobian = power_network.get_VD_jacobian()
            bus_type_init = power_network.get_bus_type_vec()
            
            print("V")
            #print(V)
            print(VD_vec)
            """

            #print(V_updated)
            #print(Q_calc)
        
        elif (i==5):
            break
        else:
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, P_updated, Q_updated = iterate_fdclf(num_buses, V_updated, V_vec_1, V_vec_2, delta_updated, delta_vec, V_vec, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            """
            Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
            power_network = PQ_to_PV(bus_type_init, Q_updated, Q_max, V_updated, power_network)
            print("V_updated")
            print(V_updated)
            V = power_network.get_V_vec()
            Q_calc = power_network.get_Q_vec()
            VD_vec, VD_jacobian = power_network.get_VD_jacobian()
            PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
            delta = power_network.get_delta_calc()
            PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
            VD_vec, VD_jacobian = power_network.get_VD_jacobian()
            bus_type_init = power_network.get_bus_type_vec()
            """
            i += 1
    return P_updated, Q_updated