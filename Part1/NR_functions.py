import numpy as np
import pandas as pd
import os
import cmath
from NR_network import Network, Buses, PQ, VD


# Function to create the Y-bus matrix
# shape needs to be specified per now. Fix this?
def Ybus(file, shape):
    df_impedances = pd.read_csv(file, sep=";")
    Z_values = np.zeros((shape,shape), dtype=complex)
    Y_bus = np.zeros((shape,shape), dtype=complex)

    for x in range(df_impedances.shape[0]):
        from_line = df_impedances['From_line'][x]
        to_line = df_impedances['To_line'][x]

        # Adding the diagonal elements to the Y-bus
        Y_bus[from_line-1][from_line-1] += 1/(df_impedances['R'][x] + df_impedances['X'][x]*1j)
        Y_bus[from_line-1][from_line-1] += 1/2*(df_impedances['Full_line_B'][x]*1j)
        Y_bus[to_line-1][to_line-1] += 1/(df_impedances['R'][x] + df_impedances['X'][x]*1j)
        Y_bus[to_line-1][to_line-1] += 1/2*(df_impedances['Full_line_B'][x]*1j)

        # Z values for off diagonal elements
        Z_values[from_line-1][to_line-1] = df_impedances['R'][x] + df_impedances['X'][x]*1j
        Z_values[to_line-1][from_line-1] = df_impedances['R'][x] + df_impedances['X'][x]*1j
    # Adding off diagonal elements
    for i in range(shape):
        for j in range(shape):
            if(Z_values[i][j] != 0. +0.j):
                if(i != j):
                    Y_bus[i][j] = - 1/Z_values[i][j]
            else:
                if(i != j):
                    Y_bus[i][j] = Z_values[i][j]
    return Y_bus


# Reading bus_data from file
# Adding the information to class objects
def read_buses(file):
    bus_vec = []
    df_buses = pd.read_csv(file, sep=";")
    for x in range(df_buses.shape[0]):
        V = df_buses['V'][x]
        delta = df_buses['Angle'][x]
        P = df_buses['P_gen'][x]
        Q = df_buses['Q_gen'][x]
        bus_type = df_buses['BusType'][x]
        bus_num = df_buses['BusNum'][x]
        new_bus = Buses(P, Q, V, delta, bus_num, bus_type)
        bus_vec.append(new_bus)
        
    return bus_vec


# Calculation of P_values
def P_Calc(V, YBus, BusNum, delta, P):
    P_Calc = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        if (np.isnan(P[i])):
            P_Calc[i] = np.nan
        else:
            for j in range (len(BusNum)):
                P_Calc[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.cos(delta[i]-delta[j])
                            +np.imag(YBus)[i][j]*cmath.sin(delta[i]-delta[j]))) 
    return P_Calc


#Function to caluculate the known Q values
def Q_Calc(V, YBus, BusNum, delta, Q):
    Q_Calc = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        if (np.isnan(Q[i])):
            Q_Calc[i] = np.nan
        else:
           for j in range (len(BusNum)):
                Q_Calc[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.sin(delta[i]-delta[j])
                                         -np.imag(YBus)[i][j]*cmath.cos(delta[i]-delta[j]))) 
    return Q_Calc


def get_PQ_calc(P_calculated, Q_calculated):
    PQ_calc = []
    for x in range(len(P_calculated)):
        if (np.isnan(P_calculated[x]) == False):
            PQ_calc.append(P_calculated[x])
            #PQ_calc[count] = P_calculated[x]

    for x in range(len(Q_calculated)):
        if (np.isnan(Q_calculated[x]) == False):
            PQ_calc.append(Q_calculated[x])
            #PQ_calc[count] = Q_calculated[x]
    return PQ_calc



def make_jacobian(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, Ybus):
    j = np.zeros((7,7), dtype=complex)
    
    for x in range(len(PQ_vec)):
        for y in range(len(PQ_vec)):
            if (PQ_jacobian[x].Bus_num == VD_jacobian[y].Bus_num):
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'D'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 0
                        else:                            
                            j[x,y] += V[PQ_jacobian[x].Bus_num]*(-np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])+np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))      
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'V'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 2*V[PQ_jacobian[x].Bus_num]*np.real(Ybus[PQ_jacobian[x].Bus_num,i])
                        else:
                            j[x,y] += V[PQ_jacobian[x].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i])+np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i]))      

                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*V[y]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))+np.imag(Ybus[PQ_jacobian[x].Bus_num,i]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])))

                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += 2*V[PQ_jacobian[x].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])-np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))
           
            else:
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*V[VD_jacobian[y].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))
                
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])+np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])))
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*V[VD_jacobian[y].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])))
    return j



def delta_VD(PQ_vec, PQ_calc, j_inv):
    delta_PQ = np.array(PQ_vec) - np.array(PQ_calc)
    delta_VD = np.matmul(j_inv,delta_PQ)
    return delta_VD




def updateVD(VD_vec, delta_vd):
    VD_vec_updated = VD_vec.copy()
    VD_vec_updated = np.array(VD_vec) + np.array(delta_vd)
    return VD_vec_updated



def updateVD_vec(VD_vec_current,delta_current,V_current):
    delta_updated = delta_current.copy()
    delta_updated = np.array(delta_updated, complex)
    V_updated = V_current.copy()
    V_updated = np.array(V_updated, complex)
    c = 0
    for x in range(len(delta_updated)):
        if (np.isnan(delta_updated[x])):
            delta_updated[x] = VD_vec_current[c]
            c += 1
    for x in range(len(V_updated)):
        if (np.isnan(V_updated[x])):
            V_updated[x] = VD_vec_current[c]
            c += 1
    return delta_updated, V_updated



def updatePQ_vec(PQ_vec, V_current, delta_current, Ybus, bus_num_init, P_init, Q_init):
    PQ_vec_updated = PQ_vec.copy() #Hvorfor gjør vi dette?
    P_calc_updated = P_Calc(V_current, Ybus, bus_num_init, delta_current, P_init)
    Q_calc_updated = Q_Calc(V_current, Ybus, bus_num_init, delta_current, Q_init)
    PQ_vec_updated = get_PQ_calc(P_calc_updated, Q_calc_updated) 
    return PQ_vec_updated