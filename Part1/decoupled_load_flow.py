import numpy as np
import pandas as pd
import cmath

from NR_functions import read_buses, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV
from DCLF_functions import Ybus_dclf, printing_Y_bus, iterate_dclf
from NR_network import Network, Buses, PQ, VD



bus_vec = read_buses('Part1/Busdata.csv')
power_network = Network(bus_vec)
#print(power_network.get_bus(2))

Ybus = Ybus_dclf('Part1/impedances.csv', 5)
convergence = 0.00001
Q_max = [0.5, 5, 1.5,5,5]


def DCLF(Ybus, power_network, convergence, Q_max):
    V_init = power_network.get_V_calc()
    delta_init = power_network.get_delta_calc()
    P_init = power_network.get_P_vec()
    Q_init = power_network.get_Q_vec()
    bus_num_init = power_network.get_bus_num_vec()
    V = power_network.get_V_vec()
    delta = power_network.get_delta_vec()
    PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
    VD_vec, VD_jacobian = power_network.get_VD_jacobian()
    bus_type_init = power_network.get_bus_type_vec()
    num_buses = len(bus_num_init)
    

    delta_vd = [1,1,1,1,1,1,1]
    PQ_vec_updated = PQ_vec.copy()
    i= 0
    #print("V")
    #print(V)
    while(abs(max(np.real(delta_vd))) > convergence):
        if (i==0):
            print("Iteration", i+1, ": \n")
            bus_type = power_network.get_bus_type_vec()
            PQ_vec_updated, delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, delta_vd = iterate_dclf(VD_jacobian, PQ_jacobian, PQ_vec, PQ_vec, num_buses, V, delta, V_init, delta_init, Ybus, bus_num_init, P_init, Q_init, VD_vec, power_network, bus_type_init, bus_type)
            i += 1
            #print("V_updated")
            #print(V_updated)
            #print("delta_updated")
            #print(delta_updated)
            P_updated = P_Updated(V_updated, Ybus, bus_num_init, delta_updated)
            Q_updated = Q_Updated(V_updated, Ybus, bus_num_init, delta_updated)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_init)
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
            PQ_vec_updated, delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, delta_vd = iterate_dclf(VD_jacobian, PQ_jacobian, PQ_vec, PQ_vec_updated, num_buses, V, delta, V_updated, delta_updated, Ybus, bus_num_init, P_calc, Q_calc, VD_vec_current,power_network, bus_type_init, bus_type)
            P_updated = P_Updated(V_updated, Ybus, bus_num_init, delta_updated)
            Q_updated = Q_Updated(V_updated, Ybus, bus_num_init, delta_updated)
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_init)
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
    return PQ_vec_updated, P_updated, Q_updated 
      


a, b, c = DCLF(Ybus, power_network, convergence, Q_max)

print(a)
print(b)
print(c)


