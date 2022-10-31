from operator import matmul
from re import I
import numpy as np
import pandas as pd
import cmath

from NR_functions import read_buses, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV
from DCLF_functions import Ybus_dclf, printing_Y_bus, iterate_dclf
from NR_network import Network, Buses, PQ, VD






def Ybus_fdclf(file, shape, bus_file, power_network):

    bus_vec = read_buses(bus_file)


    df_impedances = pd.read_csv(file, sep=";")
    Z_values = np.zeros((shape,shape), dtype=complex)
    Y_bus = np.zeros((shape,shape), dtype=complex)
    
    num_PV = 0
    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type==1):
            num_PV +=1


    for x in range(df_impedances.shape[0]):

        bus_type = bus_vec[x].bus_type
        from_line = df_impedances['From_line'][x]
        to_line = df_impedances['To_line'][x]

        # Adding the diagonal elements to the Y-bus
        Y_bus[from_line-1][from_line-1] += 1/(df_impedances['X'][x]*1j)
        Y_bus[to_line-1][to_line-1] += 1/(df_impedances['X'][x]*1j)

        # Z values for off diagonal elements
        Z_values[from_line-1][to_line-1] = df_impedances['X'][x]*1j
        Z_values[to_line-1][from_line-1] = df_impedances['X'][x]*1j
    # Adding off diagonal elements
    for i in range(shape):
        for j in range(shape):
            if(Z_values[i][j] != 0. +0.j):
                if(i != j):
                    Y_bus[i][j] = - 1/Z_values[i][j]
            else:
                if(i != j):
                    Y_bus[i][j] = Z_values[i][j]

    b_dash = Y_bus.copy()
    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0):
            b_dash = np.delete(b_dash, x, 0)
            b_dash = np.delete(b_dash, x, 1)
    
    b_dash = b_dash*-(1.j)    
    
    b_double_dash = Y_bus.copy()
    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0 or bus_vec[x].bus_type == 1):
            b_double_dash = np.delete(b_double_dash, x, 0)
            b_double_dash = np.delete(b_double_dash, x, 1)
    
    b_double_dash = b_double_dash*-(1.j) 

    print(b_dash)
    print(b_double_dash)

    return b_dash, b_double_dash



def iterate_fdclf(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q):
    
    #1 Get delta_P
    P_updated = P_Calc(V, Ybus, bus_num_init, delta, P)

    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            P_updated = np.delete(P_updated, x, 0)
 
    delta_P = P_vec_FD - P_updated
   
    #2 Get delta_Delta
    b_dash_inv =  np.linalg.inv(b_dash)
   
    delta_Delta = np.matmul(-b_dash_inv,(delta_P/V_vec_1))

    print('delta_Delta')
    print(delta_Delta)
    print('----')
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
    Q_updated  = Q_Calc(V, Ybus, bus_num_init, delta_updated, Q)

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
    print(V)
    print(delta_updated)

    return V, delta_updated, P_updated, Q_updated, V_vec_1_updated, V_vec_2_updated



    

