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



Ybus = Ybus('Part1/impedances.csv', 5)




print(Ybus)





