from operator import matmul
from re import I
import numpy as np
import pandas as pd
import cmath
from pandas import *

from Newton_raphson.NR_functions import read_buses, P_Updated, Q_Updated, Q_max_violation

def Ybus_fdclf(file, shape, bus_file, power_network, Ybus):

    bus_vec = read_buses(bus_file)
    #Making B dash matrix
    b_dash = np.zeros((shape, shape), dtype=float)
    for i in range(shape):
        for j in range(shape):
            b_dash[i][j] = np.real(np.imag(Ybus[i][j]))

    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0):
            b_dash = np.delete(b_dash, x, 0)
            b_dash = np.delete(b_dash, x, 1) 
    
    #Making B double dash matrix
    b_double_dash = np.zeros((shape, shape), dtype=float)
    for i in range(shape):
        for j in range(shape):
            b_double_dash[i][j] = np.real(np.imag(Ybus[i][j]))

    i = 0
    for x in range(len(bus_vec)):
        if(bus_vec[x].bus_type == 0 or bus_vec[x].bus_type == 1):
            b_double_dash = np.delete(b_double_dash, x-i, 0)
            b_double_dash = np.delete(b_double_dash, x-i, 1)
            i += 1    

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




def iterate_fdclf_2(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, b_dash, b_double_dash, P, Q, Q_max, power_network):
    
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
    """
    bus_type = bus_type_vec
    if(Q_violated(Q_max, Q_updated_return, bus_type)):
        Q_updated, power_network = Q_max_violation(Q_updated_return, Q_max, bus_num_init, V, power_network)
        bus_type = power_network.get_bus_type_vec()
        V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
        Q_vec_FD = power_network.get_Q_vec_FD()
        P_vec_FD = power_network.get_P_vec_FD()
        print("V_vec")
        print(V_vec_1)
        print(V_vec_2)
        
        #VD_vec_current = VD_vec_Qmax(VD_vec, VD_vec_current, bus_type, bus_type_init, V)
        
        
        V  = power_network.get_V_vec()
        #VD_vec_current = updateVD_vec(VD_vec_current, delta_vd, bus_type_init, bus_type, delta, V)


        #delta_updated, V_updated = updateVD(VD_vec_current,delta, V, bus_type_init, bus_type)
        #Q_calc = Q_calc_violated(bus_type_init,bus_type, Q_updated, Q_calc)
    """ 
    return V, delta_updated, delta_P, delta_Delta, p_updated_return, Q_updated_return, V_vec_1_updated, V_vec_2_updated, power_network
    
    


def Q_violated(Q_max, Q_Calc, bustype):
    for x in range(len(Q_max)):
        if (abs(Q_Calc[x]) > abs(Q_max[x]) and bustype[x] == 1):
            return True
    return False

#New print-functions
#Better output B_dash
def printing_B_dash(B_dash):
    df = DataFrame(B_dash)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('B dash: \n', df, '\n')
    return 
#Better output B_dobbel_dash
def printing_B_double_dash(B_double_dash):
    df = DataFrame(B_double_dash)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('B double dash: \n', df, '\n')
    return 
