from lib2to3.pgen2.pgen import DFAState
from logging import error
import numpy as np
import pandas as pd
import cmath
from pandas import *
from Newton_raphson.NR_functions import printing_Y_bus

# Function to create the Y-bus matrix
def Ybus_trans(file, shape):
    df_impedances = pd.read_csv(file, sep=";")
    Z_values = np.zeros((shape,shape), dtype=complex)
    Y_bus = np.zeros((shape,shape), dtype=complex)

    for x in range(df_impedances.shape[0]):
        from_line = df_impedances['From_line'][x]
        to_line = df_impedances['To_line'][x]

        # Adding the diagonal elements to the Y-bus
        if df_impedances['X'][x] == 0:
             Y_bus[from_line-1][from_line-1] == 0
             Y_bus[to_line-1][to_line-1] == 0
        else:
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
    
    
    read_transformers(Y_bus, "PartA/Task4/Data_transformers.csv", shape)
    
    return Y_bus

#Reading transformer alterations from file
#Alteration of Y_bus due to added transformers, both phase-shifting and tap
def read_transformers(Y_bus, file, shape):
    df_trans_info = pd.read_csv(file, sep=";")
    
    
    for x in range(df_trans_info.shape[0]):
        from_line = df_trans_info['From_line'][x]
        to_line = df_trans_info['To_line'][x]     
        X = df_trans_info['X'][x] 
        a_mag = df_trans_info['a_magnitude'][x]       
        a_angle = df_trans_info['a_angle'][x]  
        
        a = complex(df_trans_info['a_magnitude'][x], np.deg2rad(df_trans_info['a_angle'][x]))
        if df_trans_info['X'][x] == 0:
            Y_transformer = 0
        else:
            Y_transformer = 1/complex(0, df_trans_info['X'][x])
        
        Y_bus[from_line-1][from_line-1] += Y_transformer/a**2
        Y_bus[from_line-1][to_line-1] -= Y_transformer/np.conjugate(a)
        Y_bus[to_line-1][from_line-1] -= Y_transformer/a
        Y_bus[to_line-1][to_line-1] += Y_transformer
    
    return Y_bus

