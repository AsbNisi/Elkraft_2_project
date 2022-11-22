from operator import matmul
from re import I
import numpy as np
import pandas as pd
import cmath
from pandas import * 

from Newton_raphson.NR_functions import P_Updated, Q_Updated, Q_max_violation

def Ybus_fdclf(bus_vec, shape, Ybus):

    #Making B dash matrix
    b_dash = np.zeros((shape, shape), dtype=float)
    for i in range(shape):
        for j in range(shape):
            b_dash[i][j] = np.real(np.imag(Ybus[i][j]))

    for x in range(len(bus_vec)):
        if(bus_vec[x] == 0):
            b_dash = np.delete(b_dash, x, 0)
            b_dash = np.delete(b_dash, x, 1) 
    
    #Making B double dash matrix
    b_double_dash = np.zeros((shape, shape), dtype=float)
    for i in range(shape):
        for j in range(shape):
            b_double_dash[i][j] = np.real(np.imag(Ybus[i][j]))

    i = 0
    for x in range(len(bus_vec)):
        if(bus_vec[x] == 0 or bus_vec[x] == 1):
            b_double_dash = np.delete(b_double_dash, x-i, 0)
            b_double_dash = np.delete(b_double_dash, x-i, 1)
            i += 1    

    return b_dash, b_double_dash



def iterate_fdclf(num_buses, bus_num_init, V, V_vec_1, V_vec_2, delta, delta_vec, Ybus, bus_type_vec, P_vec_FD, Q_vec_FD, Q_max, power_network, method, Q_limit, reactive_limits_method):
    
    #1 Creates b_dash and b_double_dash 
    b_dash, b_double_dash = Ybus_fdclf(bus_type_vec, len(bus_type_vec), Ybus)

    #2 Get delta_P
    P_updated = P_Updated(V, Ybus, bus_num_init, delta)
    P_updated_return = P_updated.copy()
    
    #Deletes Slack-element i P_vec, to calculate delta_Delta
    P_updated = P_Updated_fixed(P_updated, bus_type_vec, num_buses)
 
    delta_P = P_vec_FD - P_updated
       
    #3 Get delta_Delta
    b_dash_inv =  np.linalg.inv(b_dash)
   
    delta_Delta = np.matmul(-b_dash_inv,(delta_P/V_vec_1))

    #Deletes Slack element in delta vec, to calculate delta updated
    delta_vec = Delta_vec_fixed(delta_vec, bus_type_vec, delta_Delta)
    
    delta_updated = delta_vec + delta_Delta
    delta_updated = delta_updated.tolist()

    #Restores delta_updated with all elements including Slack-bus element
    delta_updated = Delta_vec_restored(delta_updated, bus_type_vec, num_buses)

    #4 Find Q with new delta values
    if (method == 1):
        Q_updated  = Q_Updated(V, Ybus, bus_num_init, delta_updated)
    if (method == 2):
        Q_updated  = Q_Updated(V, Ybus, bus_num_init, delta)
    Q_updated_return = Q_updated.copy()

    #Deletes Slack-element and PV_element i Q_vec, to calculate delta_V
    Q_updated = Q_Updated_fixed(Q_updated, bus_type_vec, num_buses)

    delta_Q = Q_vec_FD - Q_updated

    #5 Find new V values
    b_double_dash_inv = np.linalg.inv(b_double_dash)

    delta_V = np.matmul(-b_double_dash_inv ,(delta_Q/V_vec_2))
    
    #Creates V_vec without Slack- and PV-elements. To calculate V_updated.
    V_limited = V_vec_fixed(V, bus_type_vec, num_buses)

    V_updated =  V_limited + delta_V 
    
    #Restores V_vec with Slack- and PV-elements. 
    V_updated = V_vec_restored(V, bus_type_vec, num_buses, V_updated)

    #Updating values for V_vec_1 and V_vec_2
    V_vec_1_updated, V_vec_2_updated = Update_V_vec(bus_type_vec, V_vec_1, V_vec_2, V_updated)
    
    #Checking for Q_max if wanted
    if (Q_limit):
        if(Q_violated(Q_max, Q_updated_return, bus_type_vec) and reactive_limits_method == 'before'):
            Q_updated, power_network = Q_max_violation(Q_updated_return, Q_max, bus_num_init, V, power_network)
            bus_type_vec = power_network.get_bus_type_vec()
            V_vec_1, V_vec_2 = power_network.get_V_vec_FD()
            Q_vec_FD = power_network.get_Q_vec_FD()
            P_vec_FD = power_network.get_P_vec_FD()
            V_vec_1_updated, V_vec_2_updated = Update_V_vec(bus_type_vec, V_vec_1, V_vec_2, V_updated)
        
     
    return V_updated, delta_updated, delta_Delta, delta_V, P_updated_return, Q_updated_return, V_vec_1_updated, V_vec_2_updated, power_network, bus_type_vec, Q_vec_FD, P_vec_FD
    
    

#Checking violation
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

#Updates V_vecs
def Update_V_vec(bus_type_vec, V_vec_1, V_vec_2, V):
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
    return V_vec_1_updated, V_vec_2_updated

#Fixing Q_updated
def Q_Updated_fixed(Q_updated, bus_type_vec, num_buses):
    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 0 or bus_type_vec[x] == 1):
            Q_updated = np.delete(Q_updated, x + i, 0)
            i -= 1
    return Q_updated

#Fixing P_updated
def P_Updated_fixed(P_updated, bus_type_vec, num_buses):
    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            P_updated = np.delete(P_updated, x, 0)
    return P_updated

#Fixing delta_vec
def Delta_vec_fixed(delta_vec, bus_type_vec, delta_Delta):
    i=0
    if(len(delta_vec) > len(delta_Delta)):
        for x in range(len(delta_vec)):
            if (bus_type_vec[x] == 0):
                    delta_vec = np.delete(delta_vec, x-i, 0)
                    i += 1 
    return delta_vec
#Restores delta_vec
def Delta_vec_restored(delta_updated, bus_type_vec, num_buses):
    for x in range(num_buses):
        if (bus_type_vec[x] == 0):
            delta_updated.insert(x,0)
    return delta_updated
    
#Fixing V_vec
def V_vec_fixed(V, bus_type_vec, num_buses):
    V_limited = []
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V_limited.append(V[x])
        else:
            continue 
    return V_limited 
#Restores V_vec
def V_vec_restored(V, bus_type_vec, num_buses, V_updated):
    i = 0
    for x in range(num_buses):
        if (bus_type_vec[x] == 2):
            V[x] = V_updated[i]
            i += 1
        else:
            continue
    return V