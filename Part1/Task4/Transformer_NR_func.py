from lib2to3.pgen2.pgen import DFAState
from logging import error
import numpy as np
import pandas as pd
import os
import cmath
from Transformer_NR_network import Network, Buses, PQ, VD
from pandas import *

#from symbol import power
def printing_Y_bus(Ybus):
    df = DataFrame(Ybus)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('Ybus: \n', df, '\n')
    return 

def printing_jacobian(j):
    df1 = DataFrame(np.real(j))
    df1.index = np.arange(1, len(df1)+1)
    df1.columns = np.arange(1, len(df1)+1)
    print('Jacobian: \n', df1, '\n')
    return 

def printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type):
    print ('Updated vales for iteration number #. Values in pu and rad.\n')
    d = {}
    for i in range (len(bus_num_init)):
        d[bus_num_init[i]+1] = bus_type[i], np.real(V_updated[i]), np.real(delta_updated[i]), np.real(P_updated[i]), np.real(Q_updated[i])
    print ("{:<7} {:<10} {:<9} {:<9} {:<14} {:<10}".format('Bus #',' Bus Type','Voltage','Angle','Active Power','Reactive Power'))    
    for k, v in d.items():
        BusType, voltage, angle, active, reactive = v
        print("{:<8} {:<9} {:<9} {:<9} {:<14} {:<16}".format(k,BusType, round(voltage,4), round(angle,4), round(active,4), round(reactive,4)))
    print('\n')
    return 


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
    printing_Y_bus(Y_bus)
    
    """
    #Alternating Y_bus due to the phase-shifting Transformer on line 1-4
    a_phase_shift = 1*(np.cos(np.deg2rad(4)) + complex(0, np.sin(np.deg2rad(4))))
    Y_phase_shift = 1/complex(0,0.2) 
    Y_bus[0][0] += Y_phase_shift/a_phase_shift**2
    Y_bus[0][3] -= Y_phase_shift/np.conjugate(a_phase_shift)
    Y_bus[3][0] -= Y_phase_shift/a_phase_shift
    Y_bus[3][3] += Y_phase_shift
    
    #Alternating Y_bus due to the phase-shifting Transformer on line 1-4
    a_tap = 0.98
    Y_tap = 1/complex(0,0.1) 
    Y_bus[3][3] += Y_tap/a_tap**2
    Y_bus[3][4] -= Y_tap/np.conjugate(a_tap)
    Y_bus[4][3] -= Y_tap/a_tap
    Y_bus[4][4] += Y_tap
    """
    printing_Y_bus(Y_bus)
    
    read_transformers(Y_bus, "Data_transformers.csv", shape)
    
    return Y_bus

#Reading transformer alterations from file
#Alteration of Y_bus due to added transformers, both phase-shifting and tap
def read_transformers(Y_bus, file, shape):
    df_trans_info = pd.read_csv(file, sep=";")
    
    for x in range(df_trans_info.shape[0]):
        from_line = df_trans_info['From_line'][x]
        to_line = df_trans_info['To_line'][x]
        
        print(to_line)
        
        print(df_trans_info['a_magnitude'][x])
        #a = complex(float(df_trans_info['a_magnitude'][x]), float(df_trans_info['a_angle'][x]))
        Y_transformer = 1/complex(0, df_trans_info['X'][x])
        
        Y_bus[from_line-1][from_line-1] += Y_transformer#/a**2
        Y_bus[from_line-1][to_line-1] -= Y_transformer#/a
        Y_bus[to_line-1][from_line-1] -= Y_transformer#/a
        Y_bus[to_line-1][to_line-1] += Y_transformer
    print("This is the new YBus")
    printing_Y_bus(Y_bus)
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
    j = np.zeros((len(PQ_vec),len(PQ_vec)), dtype=complex)
    
    for x in range(len(PQ_vec)):
        for y in range(len(PQ_vec)):
            if (PQ_jacobian[x].Bus_num == VD_jacobian[y].Bus_num):
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'D'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 0
                        else:                            
                            j[x,y] += V[PQ_jacobian[x].Bus_num]*V[i]*(-np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])+np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))      
                
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'V'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 2*V[i]*np.real(Ybus[PQ_jacobian[x].Bus_num,i])
                        else:
                            j[x,y] += V[i]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i])+np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i]))      

                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 0
                        else: 
                            j[x,y] += V[PQ_jacobian[x].Bus_num]*V[i]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))+np.imag(Ybus[PQ_jacobian[x].Bus_num,i]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])))

                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'V'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += -2*V[i]*np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i])
                        else: 
                            j[x,y] += V[i]*(np.real(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[i])-np.imag(Ybus[PQ_jacobian[x].Bus_num,i])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[i]))
            else:
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*V[VD_jacobian[y].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))
                
                if (PQ_jacobian[x].Bus_type == 'P' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*(np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])+np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num])*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])))
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*V[VD_jacobian[y].Bus_num]*(-np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num]))-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])))
    printing_jacobian(j)
    return j



def delta_VD(PQ_vec, PQ_calc, j_inv):
    delta_PQ = np.array(PQ_vec) - np.array(PQ_calc)
    delta_VD = np.matmul(j_inv,delta_PQ)
    return delta_VD




def updateVD(VD_vec, delta_vd, bus_type_init, bus_type, delta, V):
    VD_vec_updated = VD_vec.copy()
    """
    c = 0
    for x in range(len(bus_type)):
        if(np.isnan(delta[x])):
            c += 1
    for x in range(len(bus_type)):
        if(np.isnan(V[x])):
            c += 1
        if (bus_type_init[x] != bus_type[x]):
            VD_vec.insert(c-1, 1)
    print("VD_vec")
    print(VD_vec)
    """
    VD_vec_updated = np.array(VD_vec) + np.array(delta_vd)
    #VD_vec_updated = VD_vec_updated.tolist
    return VD_vec_updated



def updateVD_vec(VD_vec_current,delta,V, bus_type_init, bus_type):
    delta_current = delta.copy()
    V_current = V.copy()
    c = 0
    for x in range(len(delta_current)):
        if (np.isnan(delta_current[x])):
            delta_current[x] = VD_vec_current[c]
            c += 1

    for x in range(len(V_current)):

        #if (bus_type_init[x] != bus_type[x]):
            #V_current[x] = 1
            #VD_vec_current.insert(c,1)

        if (np.isnan(V_current[x])): #elif
            V_current[x] = VD_vec_current[c]
            c += 1     
    
    return delta_current, V_current



def updatePQ_vec(PQ_vec, V_current, delta_current, Ybus, bus_num_init, P_init, Q_init):
    #PQ_vec_updated = PQ_vec.copy() #Hvorfor gjør vi dette?
    P_calc_updated = P_Calc(V_current, Ybus, bus_num_init, delta_current, P_init)
    Q_calc_updated = Q_Calc(V_current, Ybus, bus_num_init, delta_current, Q_init)
    return get_PQ_calc(P_calc_updated, Q_calc_updated) 
    

#Function to caluculate updated P values
def P_Updated(V, YBus, BusNum, delta):
    P_Updated = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        for j in range (len(BusNum)):
            P_Updated[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.cos(delta[i]-delta[j])
                        +np.imag(YBus)[i][j]*cmath.sin(delta[i]-delta[j]))) 
    return P_Updated

#Function to caluculate updated Q values
def Q_Updated(V, YBus, BusNum, delta):
    Q_updated = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        for j in range (len(BusNum)):
            Q_updated[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.sin(delta[i]-delta[j])
                         -np.imag(YBus)[i][j]*cmath.cos(delta[i]-delta[j]))) 
    return Q_updated


def Q_max_violation(Q_updated, Q_max, bus_num, V, power_network):
    V_updated = V.copy()
    Buses = power_network.buses
    #Q_updated = [0.75, -2.37, 1.07, 2.06, -0.43]
    for i in range (len(Q_max)):
        if Q_max[i] == '':
            continue
        if Q_max[i] < Q_updated[i]:
            #print('Q_max is violated for bus ', bus_num[i+1], 'and needs to be type switched.')
            Buses[i].bus_type = 2
            #powebus_type[i] = 2
            Q_updated[i] = Q_max[i]
            V_updated[i] = np.nan
            Buses[i].Q = Q_max[i]
            Buses[i].V = np.nan
            #print(V)
            #print(bus_type)  
        else:
            #print('Bus', bus_num[i+1], 'is within its boundaries.')
            continue
        #print(power_network.get_Q_vec())
        #print(power_network.get_V_vec())
        power_network = Network(Buses)
    return Q_updated, power_network

def PQ_to_PV(bus_type_init, Q_updated, Q_max, V_updated, power_network):
    Buses = power_network.buses
    #Q_updated = [0.75, -2.37, 1.07, 2.06, -0.43]
    for i in range (len(Q_max)):
        if (bus_type_init[i] != power_network.get_bus_type_vec()[i] and Q_updated[i] < Q_max[i]):
            Buses[i].bus_type = 1
            Buses[i].Q = np.nan
            Buses[i].V = V_updated[i]
        else:
            continue
        power_network = Network(Buses)
    return power_network

def iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, PQ_vec_updated, num_buses, V, delta, V_vec, delta_vec, Ybus, bus_num_init, P_init, Q_init, VD_vec_current, power_network, bus_type_init, bus_type):
    #1 Updates V_values and delta_values in separate vectores tougether with given values.  

    delta_updated, V_updated = updateVD_vec(VD_vec_current, delta, V, bus_type_init, bus_type)

    #2 Calculates new values for P and Q separately

    P_calc = P_Calc(V_updated, Ybus, bus_num_init, delta_updated, P_init)

    Q_calc = Q_Calc(V_updated, Ybus, bus_num_init, delta_updated, Q_init)
    
    #3 Updates RHS of inverted-jacobi-equation 
    PQ_calc_updated = get_PQ_calc(P_calc, Q_calc) 

    #4 Creating new jacobian
    j = make_jacobian(VD_jacobian, PQ_jacobian, PQ_calc_updated, num_buses, V_updated, delta_updated, Ybus)
    
    #5 Inverting jacobian 
    j_inv = np.linalg.inv(j)

    #6 Updates LHS of  inverted-jacobi-equation
    delta_vd = delta_VD(PQ_vec, PQ_calc_updated, j_inv)
    
    #7 Updates values for V and delta
    VD_vec_current = updateVD(VD_vec_current, delta_vd, bus_type_init, bus_type, delta, V)
    

    #8 Updates V_values and delta_values in separate vectores tougether with given values.  
    delta_updated, V_updated = updateVD_vec(VD_vec_current,delta,V, bus_type_init, bus_type)

    #9 Updates P and Q values in one given vector 
    PQ_vec_updated = updatePQ_vec(PQ_vec, V_updated, delta_updated, Ybus, bus_num_init, P_init, Q_init)

    return PQ_vec_updated, delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, delta_vd

