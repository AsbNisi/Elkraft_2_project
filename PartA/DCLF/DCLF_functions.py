import numpy as np
import pandas as pd
import cmath 

from Newton_raphson.NR_functions import Q_violated, printing_jacobian, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, P_Updated, Q_Updated, Q_max_violation, VD_vec_Qmax, Q_calc_violated
from Newton_raphson.NR_network import Network

# Spezialized function for creating Ybus for the DCLF method

def Ybus_dclf(file, shape):
    df_impedances = pd.read_csv(file, sep=";")
    Z_values = np.zeros((shape,shape), dtype=complex)
    Y_bus = np.zeros((shape,shape), dtype=complex)

    for x in range(df_impedances.shape[0]):
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
    return Y_bus


# Function to print the Y bus
def printing_Y_bus(Ybus):
    df = pd.DataFrame(Ybus)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('Ybus: \n', df, '\n')
    return 


# Function creating the jacobian matrix for the DCLF method. Returning the nxn matrix
def make_jacobian_dclf(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, Ybus):
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
                            j[x,y] += 0
                        else:
                            j[x,y] += 0

                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    for i in range(num_buses):
                        if (i==PQ_jacobian[x].Bus_num):
                            j[x,y] += 0
                        else: 
                            j[x,y] += 0

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
                    j[x,y] += 0
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'V'):
                    j[x,y] += V[PQ_jacobian[x].Bus_num]*np.real(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.sin(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])-np.imag(Ybus[PQ_jacobian[x].Bus_num,VD_jacobian[y].Bus_num]*cmath.cos(delta[PQ_jacobian[x].Bus_num]-delta[VD_jacobian[y].Bus_num])))
                
                if (PQ_jacobian[x].Bus_type == 'Q' and VD_jacobian[y].Bus_type == 'D'):
                    j[x,y] += 0
    printing_jacobian(j)
    return j
 

# Function with all the steps for one iteration of DCLF. To be used in a loop in the main function 
def iterate_dclf(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_vec, delta_vec, Ybus, bus_num_init, P_init, Q_init, VD_vec_current, power_network, bus_type_init, Q_max, Q_limit, reactive_limits_method):
     
    #1 Calculates new values for P and Q separately
    P_calc = P_Calc(V_vec, Ybus, bus_num_init, delta_vec, P_init)
    Q_calc = Q_Calc(V_vec, Ybus, bus_num_init, delta_vec, Q_init)
    
    #2 Updates RHS of inverted-jacobi-equation 
    PQ_calc_updated = get_PQ_calc(P_calc, Q_calc) 

    #3 Creating new jacobian
    j = make_jacobian_dclf(VD_jacobian, PQ_jacobian, PQ_calc_updated, num_buses, V_vec, delta_vec, Ybus)
    
    #4 Inverting jacobian 
    j_inv = np.linalg.inv(j)

    #5 Updates LHS of  inverted-jacobi-equation
    delta_vd = delta_VD(PQ_vec, PQ_calc_updated, j_inv)
    
    #6 Updates values for V and delta
    VD_vec_current = updateVD_vec(VD_vec_current, delta_vd, bus_type_init, bus_type_init, delta, V)

    #7 Updates V_values and delta_values in separate vectores tougether with given values.  
    delta_updated, V_updated = updateVD(VD_vec_current,delta, V , bus_type_init, bus_type_init)


    #8 New P and Q vectors without nan.
    P_updated = P_Updated(V_updated, Ybus, bus_num_init, delta_updated)     #Returns P vector updated with calculated values for unknown P's instead of nan.
    Q_updated = Q_Updated(V_updated, Ybus, bus_num_init, delta_updated)     #Returns Q vector updated with calculated values for unknown Q's instead of nan.
    
    #8 Checking Q_max 
    bus_type = bus_type_init
    if (Q_limit):
        if(Q_violated(Q_max, Q_updated, bus_type) and reactive_limits_method == 'before'):
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


#PQ_to_PV
def PQ_to_PV(bus_type_init, bus_type, power_network, V_updated):
    Buses = power_network.buses
    for i in range (len(bus_type_init)):
        if (bus_type_init[i] != bus_type[i]):
            Buses[i].bus_type = 1
            Buses[i].Q = np.nan
            Buses[i].V = V_updated[i]
        else:
            continue
        power_network = Network(Buses)
    return power_network