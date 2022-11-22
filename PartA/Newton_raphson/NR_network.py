import numpy as np
import pandas as pd

# Definition of network class. A vector of buses
class Network:
    def __init__(self, buses):
        self.buses = buses

    # Class function returning vector of bus objects
    def get_buses(self):
        return self.buses()
    
    # Class function returning vector with voltage values for each bus
    def get_V_vec(self):
        V_vec = []
        for x in range(len(self.buses)):
            V_vec.append(self.buses[x].V)
        return V_vec

    # Class function returning vector with active power values for each bus
    def get_P_vec(self):
        P_vec = []
        for x in range(len(self.buses)):
            P_vec.append(self.buses[x].P)
        return P_vec
    
    # Class function returning vector with active power values for each bus
    def get_P_vec_FD(self):
        P_vec = []
        for x in range(len(self.buses)):
            if (np.isnan(self.buses[x].P) == False):
                P_vec.append(self.buses[x].P)
            else:
                continue
        return P_vec

    # Class function returning vector with reactive power values for each bus. 
    # Spezialized for Fast decoupled method
    def get_Q_vec_FD(self):
        Q_vec = []
        for x in range(len(self.buses)):
            if (self.buses[x].bus_type == 2):
                Q_vec.append(self.buses[x].Q)
            else:
                continue
        return Q_vec

    # Class function returning vector with active power values for each bus. 
    # Spezialized for Fast decoupled method
    def get_V_vec_FD(self):
        V_vec_1 = []
        V_vec_2 = []
        for x in range(len(self.buses)):
            if (self.buses[x].bus_type != 0): 
                if (self.buses[x].bus_type != 1): 
                    V_vec_1.append(1)
                else: 
                    V_vec_1.append(self.buses[x].V)
            if (self.buses[x].bus_type == 2):
                V_vec_2.append(1)
            else:
                continue
        return V_vec_1, V_vec_2

    # Class function returning vector with angle values for each bus. 
    # Spezialized for Fast decoupled method
    def get_delta_vec_FD(self):
        delta_vec = []
        for x in range(len(self.buses)):
            if (self.buses[x].bus_type != 0):
                if (np.isnan(self.buses[x].delta) == True):
                    delta_vec.append(0)
                if (np.isnan(self.buses[x].delta) == False):
                    delta_vec.append(self.buses[x].delta)
            else:
                continue
        return delta_vec

    # Class function returning vector with reactive power values for each bus. 
    def get_Q_vec(self):
        Q_vec = []
        for x in range(len(self.buses)):
            Q_vec.append(self.buses[x].Q)
        return Q_vec

    # Class function returning vector with angle values for each bus. 
    def get_delta_vec(self):
        delta_vec = []
        for x in range(len(self.buses)):
            delta_vec.append(self.buses[x].delta)
        return delta_vec

    # Class function returning vector with bus types for each bus. 
    def get_bus_type_vec(self):
        bus_type_vec = []
        for x in range(len(self.buses)):
            bus_type_vec.append(self.buses[x].bus_type)
        return bus_type_vec

    # Class function returning vector with bus number for each bus. 
    def get_bus_num_vec(self):
        bus_num_vec = []
        for x in range(len(self.buses)):
            bus_num_vec.append(self.buses[x].bus_num)
        return bus_num_vec
    
    # Class function returning vector with voltage values for first iteration. 
    # If the value is NaN, an initial guess of 1 pu is added
    def get_V_calc(self):
        V_vec = []
        for x in range(len(self.buses)):
            if(np.isnan(self.buses[x].V)):
                V_vec.append(1)
            else:
                V_vec.append(self.buses[x].V)
        return V_vec

    # Class function returning vector with delta values for first iteration. 
    # If the value is NaN, an initial guess of 0 is added
    def get_delta_calc(self):
        delta_vec = []
        for x in range(len(self.buses)):
            if (np.isnan(self.buses[x].delta)):
                delta_vec.append(0)
            else:
                delta_vec.append(self.buses[x].delta)
        return delta_vec
    
    # Class function returning vector with the known active and reactive power values 
    def get_PQ_vec(self):
        PQ_vec = []
        PQ_jacobian = []
        for x in range(len(self.buses)):
            if (np.isnan(self.buses[x].P) == False):
                PQ_vec.append(self.buses[x].P)
                businfo = PQ("P", self.buses[x].bus_num)
                PQ_jacobian.append(businfo)
        for x in range(len(self.buses)):
            if (np.isnan(self.buses[x].Q) == False):
                PQ_vec.append(self.buses[x].Q)
                businfo = PQ("Q", self.buses[x].bus_num)
                PQ_jacobian.append(businfo)
        return PQ_vec, PQ_jacobian

    # Class function returning vector with voltages and delta to be used for the calculation of the jacobian
    def get_VD_jacobian(self):
        VD_jacobian = []
        VD_vec = []

        for x in range(len(self.buses)):
            if (np.isnan(self.buses[x].delta)):
                VD_vec.append(0)
                businfo = VD("D", self.buses[x].bus_num)
                VD_jacobian.append(businfo)  

        for x in range(len(self.buses)):    
            if (np.isnan(self.buses[x].V)):
                VD_vec.append(1)
                businfo = VD("V", self.buses[x].bus_num)
                VD_jacobian.append(businfo)

        return VD_vec, VD_jacobian

# Bus class storing values for each single bus. NaN values if unknown
class Buses:
    def __init__(self, P, Q, V, delta, bus_num, bus_type):
        self.P = P
        self.Q = Q
        self.V = V
        self.delta = delta
        self.bus_num = bus_num 
        self.bus_type = bus_type

# Class helping to keep track of which values we should take the derivative of in the Jacobian
class PQ:
    def __init__(self, Bus_type, Bus_num):
        self.Bus_type = Bus_type
        self.Bus_num = Bus_num
    def get_bus_num(self):
        return self.Bus_num

# Class helping to keep track of which values we should take the derivative of in the Jacobian
class VD:
    def __init__(self, Bus_type, Bus_num):
        self.Bus_type = Bus_type
        self.Bus_num = Bus_num





