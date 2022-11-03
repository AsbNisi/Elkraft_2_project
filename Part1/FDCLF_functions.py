from operator import matmul
from re import I
import numpy as np
import pandas as pd
import cmath

from NR_functions import read_buses, Ybus, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV
from DCLF_functions import Ybus_dclf, printing_Y_bus, iterate_dclf
from NR_network import Network, Buses, PQ, VD


def Ybus_fdclf(file, shape, bus_file, power_network, Ybus):

    bus_vec = read_buses(bus_file)

    
    b_dash = Ybus.copy()
    
    for i in range(shape):
        for j in range(shape):
            b_dash[i][j] = np.imag(b_dash[i][j])

    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0):
            b_dash = np.delete(b_dash, x, 0)
            b_dash = np.delete(b_dash, x, 1)
    
    #b_dash = b_dash*-(1.j)    
    
    b_double_dash = Ybus.copy()
    for i in range(shape):
        for j in range(shape):
            b_double_dash[i][j] = np.imag(b_double_dash[i][j])

    i = 0
    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0 or bus_vec[x].bus_type == 1):
            b_double_dash = np.delete(b_double_dash, x-i, 0)
            b_double_dash = np.delete(b_double_dash, x-i, 1)
            i += 1
    
    #b_double_dash = b_double_dash*-(1.j) 

    print('B_dash')
    print(b_dash)
    print(b_dash[0][0])
    print('.........')
    print('B_double_dash')
    print(b_double_dash)
    print('.........')

    return b_dash, b_double_dash


def iterate_fdclf_1(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q):
    
    #1 Get delta_P
    P_updated = P_Updated(V, Ybus, bus_num_init, delta)
    p_updated_return = P_updated.copy()

    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            P_updated = np.delete(P_updated, x, 0)
 
    delta_P = P_vec_FD - P_updated
       
    #2 Get delta_Delta
    b_dash_inv =  np.linalg.inv(b_dash)
   
    delta_Delta = np.matmul(-b_dash_inv,(delta_P/V_vec_1))
    delta_updated = delta_vec.copy()
    
    i=0
    if(len(delta_vec) > len(delta_Delta)):
        for x in range(len(delta_vec)):
            if (bus_type_vec[x] == 0):
                    delta_vec = np.delete(delta_vec, x-i, 0)
                    i += 1 

    delta_updated = delta_vec + delta_Delta
    delta_updated = delta_updated.tolist()

    
    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            delta_updated.insert(x,0)

    #delta_updated = np.array(delta_updated)
    #3 Find Q with new delta values
    Q_updated  = Q_Updated(V, Ybus, bus_num_init, delta_updated)
    Q_updated_return = Q_updated.copy()

    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 0 or bus_type_vec[x] == 1):
            Q_updated = np.delete(Q_updated, x + i, 0)
            i -= 1

    delta_Q = Q_vec_FD - Q_updated

    #4 Find new V values
    b_double_dash_inv = np.linalg.inv(b_double_dash)
    delta_V = np.matmul(-b_double_dash_inv ,(delta_Q/V_vec_2))
    V_limited = []
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V_limited.append(V[x])
        else:
            continue 

    V_updated =  V_limited + delta_V 
    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V[x] = V_updated[i]
            i += 1
        else:
            continue

    # Updating values for V_vec_1 and V_vec_2
    V_vec_1_updated = V_vec_1.copy()
    V_vec_2_updated = V_vec_2.copy()

    i = 0
    j = 0
    for x in range(len(bus_type_vec)):
        if (bus_type_vec[x] != 0): 
            V_vec_1_updated[x-i] = V[x]
        else:
            i += 1
        if (bus_type_vec[x] == 2):
            V_vec_2_updated[x-j] = V[x]
        else:
            j += 1
    return V, delta_updated, delta_P, delta_Delta, p_updated_return, Q_updated_return, V_vec_1_updated, V_vec_2_updated




def iterate_fdclf_2(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q):
    
    #1 Get delta_P
    P_updated = P_Updated(V, Ybus, bus_num_init, delta)
    p_updated_return = P_updated.copy()

    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            P_updated = np.delete(P_updated, x, 0)
 
    delta_P = P_vec_FD - P_updated
       
    #2 Get delta_Delta
    b_dash_inv =  np.linalg.inv(b_dash)
   
    delta_Delta = np.matmul(-b_dash_inv,(delta_P/V_vec_1))
    delta_updated = delta_vec.copy()
    
    i=0
    if(len(delta_vec) > len(delta_Delta)):
        for x in range(len(delta_vec)):
            if (bus_type_vec[x] == 0):
                    delta_vec = np.delete(delta_vec, x-i, 0)
                    i += 1 

    delta_updated = delta_vec + delta_Delta
    delta_updated = delta_updated.tolist()

    
    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            delta_updated.insert(x,0)

    #delta_updated = np.array(delta_updated)
    #3 Find Q with new delta values
    Q_updated  = Q_Updated(V, Ybus, bus_num_init, delta)
    Q_updated_return = Q_updated.copy()

    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 0 or bus_type_vec[x] == 1):
            Q_updated = np.delete(Q_updated, x + i, 0)
            i -= 1

    delta_Q = Q_vec_FD - Q_updated

    #4 Find new V values
    b_double_dash_inv = np.linalg.inv(b_double_dash)
    delta_V = np.matmul(-b_double_dash_inv ,(delta_Q/V_vec_2))
    V_limited = []
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V_limited.append(V[x])
        else:
            continue 

    V_updated =  V_limited + delta_V 
    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V[x] = V_updated[i]
            i += 1
        else:
            continue

    # Updating values for V_vec_1 and V_vec_2
    V_vec_1_updated = V_vec_1.copy()
    V_vec_2_updated = V_vec_2.copy()

    i = 0
    j = 0
    for x in range(len(bus_type_vec)):
        if (bus_type_vec[x] != 0): 
            V_vec_1_updated[x-i] = V[x]
        else:
            i += 1
        if (bus_type_vec[x] == 2):
            V_vec_2_updated[x-j] = V[x]
        else:
            j += 1
    return V, delta_updated, delta_P, delta_Delta, p_updated_return, Q_updated_return, V_vec_1_updated, V_vec_2_updated

    

