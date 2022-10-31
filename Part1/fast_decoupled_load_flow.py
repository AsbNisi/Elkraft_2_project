import numpy as np
import pandas as pd
import cmath

from FDCLF_functions import Ybus_fdclf, iterate_fdclf
from NR_functions import read_buses, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV, Ybus
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
    delta = power_network.get_delta_vec()
    delta = [0,0,0,0,0]
    P = power_network.get_P_vec()
    Q = power_network.get_Q_vec()



    i= 0
    while(True):
        if (i==0):
            print("Iteration", i+1, ": \n")
            print(V)
            print(Ybus)
            print(bus_num_init)
            print(delta)
            print(P)
            bus_type = power_network.get_bus_type_vec()
            V_updated, delta_updated, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated = iterate_fdclf(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec_init, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q)
                                                                             
            i += 1
            #printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
        
        elif (i==5):
            break
        else:
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated = iterate_fdclf(num_buses, bus_num_init, V_updated, V_vec_1_updated, V_vec_2_updated, delta, delta_updated, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q)
            #printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)

            print(V_updated)
            print('delta_updated')
            print(delta_updated)
            print(P_updated)
            print(Q_updated)
            i += 1
    return P_updated, Q_updated

P, Q = FDCLF(Ybus, power_network, convergence, Q_max)

print(P)
print(Q)

