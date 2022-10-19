import numpy as np
import pandas as pd
import os
import cmath



from NR_functions import Ybus, read_buses, P_calc, Q_calc, get_PQ_calc, make_jacobian, delta_VD, updateVD, updateVD_vec, updatePQ_vec 
from NR_classdef import Network, Buses, PQ, VD

bus_vec = read_buses('Busdata.csv')

power_network = Network(bus_vec)

def NR(Ybus, V, delta, power_network ):
    
    
    num_buses = len(bus_num_init)



