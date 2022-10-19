import numpy as np
import pandas as pd
import os
import cmath



from NR_functions import Ybus, read_buses, P_Calc, Q_Calc, get_PQ_calc, make_jacobian, delta_VD, updateVD, updateVD_vec, updatePQ_vec 
from NR_network import Network, Buses, PQ, VD

bus_vec = read_buses('Part1/Busdata.csv')

power_network = Network(bus_vec)

Ybus = Ybus('Part1/impedances.csv', 5)

def NR(Ybus, power_network):
    V_init = power_network.get_V_calc()
    delta_init = power_network.get_delta_calc()
    P_init = power_network.get_P_vec()
    Q_init = power_network.get_Q_vec()
    bus_num_init = power_network.get_bus_num_vec()
    V = power_network.get_V_vec()
    delta = power_network.get_delta_vec()
    PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
    VD_vec, VD_jacobian = power_network.get_VD_jacobian()

    num_buses = len(bus_num_init)

    P_calc = P_Calc(V_init, Ybus, bus_num_init, delta_init, P_init)
    Q_calc = Q_Calc(V_init, Ybus, bus_num_init, delta_init, Q_init)

    PQ_calc = get_PQ_calc(P_calc, Q_calc)  

    j = make_jacobian(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V_init, delta_init, Ybus)
    j_inv = np.linalg.inv(j)

    delta_vd = delta_VD(PQ_vec, PQ_calc, j_inv)

    VD_vec_current = updateVD(VD_vec, delta_vd)

    delta_updated, V_updated = updateVD_vec(VD_vec_current,delta,V)

    PQ_vec_updated = updatePQ_vec(PQ_vec, V_updated, delta_updated, Ybus, bus_num_init, P_init, Q_init)

    return PQ_vec_updated 
        


a = NR(Ybus, power_network)

print(a)

