from lib2to3.pgen2.pgen import DFAState
from logging import error
import numpy as np
import pandas as pd
import cmath
from Newton_raphson.NR_network import Network, Buses
from pandas import *

#Better output Y_bus
def printing_Y_bus(Ybus):
    df = DataFrame(Ybus)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('Ybus: \n', df, '\n')
    return 
#Better output jacobian 
def printing_jacobian(j):
    df1 = DataFrame(np.real(j))
    df1.index = np.arange(1, len(df1)+1)
    df1.columns = np.arange(1, len(df1)+1)
    print('Jacobian: \n', df1, '\n')
    return 
#Better output buses
def printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type):
    V_base = 132 #kV
    S_base = 100 #MW
    print("Updated bus values:", '\n')
    bus_dict = {0: "Slack", 1: "Generator", 2: "Load"}
    
    d = {}
    for i in range (len(bus_num_init)):
        d[bus_num_init[i]+1] = bus_dict.get(bus_type[i]), np.real(V_updated[i]), np.real(V_updated[i])*V_base,np.real(delta_updated[i]), np.real(P_updated[i]), np.real(Q_updated[i]), np.real(P_updated[i])*S_base, np.real(Q_updated[i])*S_base
    print ("{:<7} {:<11} {:<12} {:<12} {:<12} {:<17} {:<13} {:<17} {:<13}".format('Bus #',' Bus Type','Voltage [pu] ', 'Voltage [kV] ','Angle [rad] ','Active Power [pu] ','Reactive Power [pu] ','Active Power [MW] ','Reactive Power [MVar] '))    
    for k, v in d.items():
        BusType, voltage_pu, voltage_act, angle, active_pu, reactive_pu, active_act, reactive_act = v
        print("{:<8} {:<10} {:<13} {:<13} {:<12} {:<18} {:<20} {:<18} {:<22}".format(k,BusType, round(voltage_pu,4),round(voltage_act,4), round(angle,4), round(active_pu,4), round(reactive_pu,4), round(active_act,4), round(reactive_act,4)))
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
    return Y_bus


# Reading bus_data from file
# Adding the information to class objects Bus
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


# Calculation of P_values, returns nan if value of P is unkonwn 
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


#Function to caluculate the known Q values, returns nan if value of Q is unkonwn 
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


#Returns the caculated part of delta_PQ (LHS of jacobi equation). The given values for P and Q is obtained by a get function in the Network class. get_PQ_vec(self)
def get_PQ_calc(P_calculated, Q_calculated):
    PQ_calc = []
    for x in range(len(P_calculated)):
        if (np.isnan(P_calculated[x]) == False):
            PQ_calc.append(P_calculated[x])

    for x in range(len(Q_calculated)):
        if (np.isnan(Q_calculated[x]) == False):
            PQ_calc.append(Q_calculated[x])
    return PQ_calc


#Creates the jacobi-matrix
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


#Obtaining LHS of inverse jacobi equation 
def delta_VD(PQ_vec, PQ_calc, j_inv):
    delta_PQ = np.array(PQ_vec) - np.array(PQ_calc)
    delta_VD = np.matmul(j_inv,delta_PQ)
    return delta_VD



#Updates values for unknown Voltages and Deltas 
def updateVD_vec(VD_vec, delta_vd, bus_type_init, bus_type, delta, V):
    VD_vec_updated = VD_vec.copy()
    VD_vec = np.array(VD_vec).tolist()
    delta_vd = np.array(delta_vd).tolist()
    c = 0
    for x in range(len(bus_type)):
        if(np.isnan(delta[x])):
            c += 1
    for x in range(len(bus_type)):
        if(np.isnan(V[x])):
            c += 1
        if (bus_type_init[x] != bus_type[x]):
            delta_vd.insert(c-1, 0)
    VD_vec_updated = np.array(VD_vec) + np.array(delta_vd)
    return VD_vec_updated


#Splits the updated values of V and delta in two separate vectors  
def updateVD(VD_vec_current, delta, V, bus_type_init, bus_type):
    
    delta_current = delta.copy()
    VD_vec_current = np.array(VD_vec_current).tolist()
    V_current = V.copy()
    c = 0
    for x in range(len(delta_current)):
        if (np.isnan(delta_current[x])):
            delta_current[x] = VD_vec_current[c]
            c += 1

    for x in range(len(V_current)):

        if (bus_type_init[x] != bus_type[x]):

            #V_current[x] = VD_vec_current[c]   # test
            c +=1
            V_current[x] = 1  #Original            
            #Per nå hardkodet. Fiks hvis tid

        elif (np.isnan(V_current[x])): 
            V_current[x] = VD_vec_current[c]
            c += 1     

    return delta_current, V_current


#Updates The calculated values for P and Q with new delta and V_values 
def updatePQ_vec(PQ_vec, V_current, delta_current, Ybus, bus_num_init, P_init, Q_init):
    P_calc_updated = P_Calc(V_current, Ybus, bus_num_init, delta_current, P_init)
    Q_calc_updated = Q_Calc(V_current, Ybus, bus_num_init, delta_current, Q_init)
    return get_PQ_calc(P_calc_updated, Q_calc_updated) 
    

#Function to caluculate updated P values. Calculates all values 
def P_Updated(V, YBus, BusNum, delta):
    P_Updated = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        for j in range (len(BusNum)):
            P_Updated[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.cos(delta[i]-delta[j])
                        +np.imag(YBus)[i][j]*cmath.sin(delta[i]-delta[j]))) 
    return P_Updated

#Function to caluculate updated Q values. Calculates all values
def Q_Updated(V, YBus, BusNum, delta):
    Q_updated = np.zeros(len(BusNum), dtype=complex)
    for i in range (len(BusNum)):
        for j in range (len(BusNum)):
            Q_updated[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.sin(delta[i]-delta[j])
                         -np.imag(YBus)[i][j]*cmath.cos(delta[i]-delta[j]))) 
    return Q_updated

#Function to check if Q is violated. Updates power network
def Q_max_violation(Q_updated, Q_max, bus_num, V, power_network):
    #V_updated = V.copy()
    Buses = power_network.buses
    #Q_updated = [0.75, -2.37, 1.07, 2.06, -0.43]



    for i in range (len(Q_max)):
        if Q_max[i] == '':
            continue
        if abs(Q_max[i]) < abs(Q_updated[i]):
            print('Q_max is violated for bus ', bus_num[i+1], 'and needs to be type switched.')
            Buses[i].bus_type = 2
            Q_updated[i] = Q_max[i]
            #V_updated[i] = np.nan
            Buses[i].Q = Q_max[i]
            Buses[i].V = np.nan  
        else:
            #print('Bus', bus_num[i+1], 'is within its boundaries.')
            continue
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

def iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_vec, delta_vec, Ybus, bus_num_init, P_init, Q_init, VD_vec_current, power_network, bus_type_init, Q_max, Q_limit):
     
    #1 Calculates new values for P and Q separately
    P_calc = P_Calc(V_vec, Ybus, bus_num_init, delta_vec, P_init)
    Q_calc = Q_Calc(V_vec, Ybus, bus_num_init, delta_vec, Q_init)
    
    #2 Updates RHS of inverted-jacobi-equation 
    PQ_calc_updated = get_PQ_calc(P_calc, Q_calc) 

    #3 Creating new jacobian
    j = make_jacobian(VD_jacobian, PQ_jacobian, PQ_calc_updated, num_buses, V_vec, delta_vec, Ybus)
    
    #4 Inverting jacobian 
    j_inv = np.linalg.inv(j)

    #5 Updates LHS of  inverted-jacobi-equation
    delta_vd = delta_VD(PQ_vec, PQ_calc_updated, j_inv)
    
    #6 Updates values for V and delta
    VD_vec_current = updateVD_vec(VD_vec_current, delta_vd, bus_type_init, bus_type_init, delta, V)
    
    #print('VD_vec_current')
    #print(VD_vec_current)
    #VD_vec_current = insert_VD_vec(delta, delta_updated, V, V_updated, VD_vec_current)
    #print(VD_vec_current)

    #7 Updates V_values and delta_values in separate vectores tougether with given values.  
    delta_updated, V_updated = updateVD(VD_vec_current,delta, V , bus_type_init, bus_type_init)


    #8 New P and Q vectors without nan.
    P_updated = P_Updated(V_updated, Ybus, bus_num_init, delta_updated)     #Returns P vector updated with calculated values for unknown P's instead of nan.
    Q_updated = Q_Updated(V_updated, Ybus, bus_num_init, delta_updated)     #Returns Q vector updated with calculated values for unknown Q's instead of nan.
    
    #8 Checking Q_max 
    bus_type = bus_type_init
    if (Q_limit):
        if(Q_violated(Q_max, Q_updated, bus_type)):
            Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
            bus_type = power_network.get_bus_type_vec()
            VD_vec, VD_jacobian = power_network.get_VD_jacobian()
            PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
            
            
            VD_vec_current = VD_vec_Qmax(VD_vec, VD_vec_current, bus_type, bus_type_init, V)
            
            
            V  = power_network.get_V_vec()
            VD_vec_current = updateVD_vec(VD_vec_current, delta_vd, bus_type_init, bus_type, delta, V)


            delta_updated, V_updated = updateVD(VD_vec_current,delta, V, bus_type_init, bus_type)
            Q_calc = Q_calc_violated(bus_type_init,bus_type, Q_updated, Q_calc)

    return delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V 




def insert_VD_vec(delta, delta_updated, V, V_updated, VD_vec):
    c = 0 
    for i in range(len(delta)):
        if(np.isnan(delta[i]) == True):
            VD_vec[c] = delta_updated[i]
            c +=1 
    for j in range(len(V)):
        if(np.isnan(V[j]) == True):
            VD_vec[c] = V_updated[j]
            c+=1
    return VD_vec


def Q_violated(Q_max, Q_Calc, bustype):
    for x in range(len(Q_max)):
        if (abs(Q_Calc[x]) > abs(Q_max[x]) and bustype[x] == 1):
            return True
    return False
        

def Q_calc_violated(bus_type_init, bus_type, Q_max, Q_calc):
    for x in range(len(bus_type)):
        if (bus_type_init[x] != bus_type[x]):
            Q_calc[x] = Q_max[x]
    return Q_calc


def VD_vec_Qmax(VD_vec, VD_vec_current, bus_type, bus_type_init, V):
    i = 0
    j = 0
    for x in range(len(VD_vec)):
        if (VD_vec[x]== 0):
            VD_vec[x] = VD_vec_current[x]
        else:
            if(bus_type[j] != bus_type_init[j]):
                VD_vec[x] = V[i]
                i += 1
                j += 1
            else:
                VD_vec[x] = VD_vec_current[x-i]
                j += 1
    return VD_vec



#Printing load flow in lines
def printing_lines(bus_vec, file, V, Ybus):
    
    #Calculate complex power flow in lines
   
    S_base = 100 #MW
    
    S_ik = np.zeros((len(bus_vec), len(bus_vec)), dtype=complex)
    for i in range(len(bus_vec)):
        for k in range(len(bus_vec)):
            S_ik[i][k] = V[i]*np.conjugate(-Ybus[i][k]*(V[i]-V[k]))*S_base
            

    print("Updated line info:", '\n')
    df_lines = pd.read_csv(file, sep=";")
    d = {}
    for i in range (len(bus_vec)):
        from_line = df_lines["From_line"][i]
        to_line = df_lines["To_line"][i]

        line = str(from_line) + " - " + str(to_line)
        
        d[line] = S_ik[from_line-1,to_line-1], np.real(S_ik[from_line-1,to_line-1]), np.imag(S_ik[from_line-1,to_line-1]), np.real(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1]), np.imag(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1]) 
    print ("{:<7} {:<23} {:<12} {:<12} {:<12} {:<12}".format('Line','Complex Power Flow [MVA] ','Active Power Flow [MW] ', 'Reactive Power Flow [MVar] ', "Active Power Loss [MW] ", "Reactive Power Loss [MVar] "))    
    for k, v in d.items():
        apparent, active, reactive, ploss, qloss = v
        print("{:<7} {:<25} {:<23} {:<27} {:<23} {:<19}".format(k, round(apparent,4),round(active,4), round(reactive,4), round(ploss,4), round(qloss,4)))
    print('\n')
    return 
