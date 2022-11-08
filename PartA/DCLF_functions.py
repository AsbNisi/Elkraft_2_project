import numpy as np
import pandas as pd
import cmath

from NR_functions import read_buses, printing_jacobian, P_Calc, Q_Calc, get_PQ_calc, delta_VD, updateVD, updateVD_vec, updatePQ_vec, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV


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
    printing_Y_bus(Y_bus)
    return Y_bus



def printing_Y_bus(Ybus):
    df = pd.DataFrame(Ybus)
    df.index = np.arange(1, len(df)+1)
    df.columns = np.arange(1, len(df)+1)
    print('Ybus: \n', df, '\n')
    return 



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



def iterate_dclf(VD_jacobian, PQ_jacobian, PQ_vec, PQ_vec_updated, num_buses, V, delta, V_vec, delta_vec, Ybus, bus_num_init, P_init, Q_init, VD_vec_current, power_network, bus_type_init, bus_type):
    #1
    print("VD_vec")
    print(VD_vec_current)
    print("delta")
    print(delta)
    print("V")
    print(V)
    delta_updated, V_updated = updateVD_vec(VD_vec_current, delta, V, bus_type_init, bus_type)
    #print("delta_U")
    #print(delta_updated)
    #print("V_U")
    #print(V_updated)
    #2
    P_calc = P_Calc(V_updated, Ybus, bus_num_init, delta_updated, P_init)
    Q_calc = Q_Calc(V_updated, Ybus, bus_num_init, delta_updated, Q_init)
    
    #3
    PQ_calc_updated = get_PQ_calc(P_calc, Q_calc) 

    #4
    j = make_jacobian_dclf(VD_jacobian, PQ_jacobian, PQ_calc_updated, num_buses, V_updated, delta_updated, Ybus)
    
    #5
    j_inv = np.linalg.inv(j)

    #6
    delta_vd = delta_VD(PQ_vec, PQ_calc_updated, j_inv)
    
    #7
    print("Vd_vec")
    print(VD_vec_current)
    print("delta_vd")
    print(delta_vd)
    VD_vec_current = updateVD(VD_vec_current, delta_vd, bus_type_init, bus_type, delta, V)
    

    #8
    delta_updated, V_updated = updateVD_vec(VD_vec_current,delta,V, bus_type_init, bus_type)
    #print(delta_updated)
    #print(V_updated)

    #9
    PQ_vec_updated = updatePQ_vec(PQ_vec, V_updated, delta_updated, Ybus, bus_num_init, P_init, Q_init)

    return PQ_vec_updated, delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, delta_vd


