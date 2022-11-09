import numpy as np
import pandas as pd
import cmath

from FDCLF_functions import Ybus_fdclf, iterate_fdclf_1, iterate_fdclf_2
from Newton_raphson.NR_functions import read_buses, Ybus
from Newton_raphson.NR_network import Network


bus_vec = read_buses('PartA/Busdata.csv')
power_network = Network(bus_vec)
Ybus = Ybus('PartA/impedances.csv', 5)

convergence = 0.00001
Q_max = [0.5, 5, 1.5,5,5]


def FDCLF(Ybus, power_network, convergence, Q_max):
    V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
    Q_vec_FD = power_network.get_Q_vec_FD()
    P_vec_FD = power_network.get_P_vec_FD()
    b_dash, b_double_dash = Ybus_fdclf('PartA/impedances.csv', 5, 'PartA/Busdata.csv', power_network, Ybus)
    bus_num_init = power_network.get_bus_num_vec()
    num_buses = len(bus_num_init)
    delta_vec_init = power_network.get_delta_vec_FD()
    bus_type_vec = power_network.get_bus_type_vec()
    V = power_network.get_V_calc()
    delta = power_network.get_delta_vec()
    delta = np.zeros(len(bus_vec))
    P = power_network.get_P_vec()
    Q = power_network.get_Q_vec()
    print(Ybus)

    delta_Delta = [1]
    delta_P = [1]
    i= 0
    while((abs(max(np.real(delta_Delta))) > convergence) and (abs(max(np.real(delta_P))) > convergence)):
        if (i==0):
            print("Iteration", i+1, ": \n")
            print('V')
            print(V)
            print('.........')
            print('delta')
            print(delta)
            print('.........')
            print('P')
            print(P)
            print('.........')
            V_updated, delta_updated, delta_P, delta_Delta, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network = iterate_fdclf_2(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec_init, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q, Q_max, power_network)
                                                                             
            i += 1
            #printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
        
        elif (i==20):
            break
        else:
            print("Iteration", i+1, ": \n")
            V_updated, delta_updated, delta_P, delta_Delta, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated, power_network = iterate_fdclf_2(num_buses, bus_num_init, V_updated, V_vec_1_updated, V_vec_2_updated, delta_updated, delta_updated, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q, Q_max, power_network)
            #printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_vec)
            
            print('V_updated')
            print(V_updated)
            print('.........')
            print('delta_updated')
            print(delta_updated)
            print('.........')
            print('P_updated')
            print(P_updated)
            print('.........')
            print('Q_updated')
            print(Q_updated)
            print('.........')
            i += 1
    return P_updated, Q_updated

P, Q = FDCLF(Ybus, power_network, convergence, Q_max)


