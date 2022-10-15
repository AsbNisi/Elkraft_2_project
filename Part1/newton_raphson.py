import numpy as np
import pandas as pd
import os
import cmath


# Reading bus data intp pandas 


def buses(file):
    df_buses = pd.read_csv(file, sep=";")
    V = df_buses.V.to_numpy()
    delta = df_buses.Angle.to_numpy()
    P = df_buses.P_gen.to_numpy()
    Q = df_buses.Q_gen.to_numpy()
    return V,delta,P,Q


V,delta,P,Q = buses('Part1/Busdata.csv')

print(V)
print(delta)
print(P)
print(Q)




# Creating Y-bus matrix

def Ybus(file, shape):
    df_impedances = pd.read_csv(file, sep=";")
    print(df_impedances.head())
    Z_values = np.zeros((shape,shape), dtype=complex)
    Y_bus = np.zeros((shape,shape), dtype=complex)

    for x in range(df_impedances.shape[0]):
        from_line = df_impedances['From_line'][x]
        to_line = df_impedances['To_line'][x]
        Z_values[from_line-1][to_line-1] = df_impedances['R'][x] + df_impedances['X'][x]*1j

    
    for i in range(shape):
        for j in range(shape):
            if(Z_values[i][j] != 0. +0.j):
                Y_bus[i][j] = 1/Z_values[i][j]
            else:
                Y_bus[i][j] = Z_values[i][j]

    for x in range(df_impedances.shape[0]):
        from_line = df_impedances['From_line'][x]
        to_line = df_impedances['To_line'][x]
        
        Y_bus[from_line-1][from_line-1] += - 1/(df_impedances['R'][x] + df_impedances['X'][x]*1j)
        Y_bus[from_line-1][from_line-1] += 1/(2*(df_impedances['Full_line_B'][x]*1j))
        Y_bus[to_line-1][to_line-1] += 1/(2*(df_impedances['Full_line_B'][x]*1j))

    return Y_bus

def P_Calc(V, YBus, BusNum, delta, P):
    P_Calc = np.zeros(BusNum, dtype=float)
    for i in range (BusNum-1):
        if P[i] == 0:
            P_Calc[i] == 0
        else:
            for j in range (BusNum-1):
                P_Calc[i] = P_Calc[i] + (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.cos(delta[i]-delta[j])
                            +np.imag(YBus)[i][j]*cmath.sin(delta[i]-delta[j]))) 
    return P_Calc

def Q_Calc(V, YBus, BusNum, delta, Q):
    Q_Calc = np.zeros(BusNum, dtype=float)
    for i in range (BusNum-1):
        if Q[i] == 0:
            Q_Calc[i] == 0
        else:
           for j in range (BusNum-1):
                Q_Calc[i] += (V[i]*V[j]*(np.real(YBus)[i][j]*cmath.sin(delta[i]-delta[j])
                                         -np.imag(YBus)[i][j]*cmath.cos(delta[i]-delta[j]))) 
    return Q_Calc



P_calc = (P_Calc(V, Ybus, 5, delta, P))
Q_calc = (Q_Calc(V, Ybus, 5, delta, Q))

print(P_calc)
print(Q_calc)

Ybus = Ybus('Part1/impedances.csv', 5)




print(Ybus)





