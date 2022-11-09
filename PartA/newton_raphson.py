import numpy as np


from NR_functions import Ybus, read_buses, insert_VD_vec, iterate_NR, P_Updated, Q_Updated, Q_max_violation, printing_buses, PQ_to_PV
from NR_network import Network

bus_vec = read_buses('PartA/Busdata.csv')

power_network = Network(bus_vec)

Ybus = Ybus('PartA/impedances.csv', 5)
convergence = 0.00001
Q_max = [0.5,5,-1.5,5,5]

def NR(Ybus, power_network, convergence, Q_max):
    V_init = power_network.get_V_calc()   #Appends 1 if nan. Otherwise given value 
    delta_init = power_network.get_delta_calc()  #Appends 0 if nan. Otherwise given value 
    P_init = power_network.get_P_vec() #Returns P_values. If unknown, returns nan
    Q_init = power_network.get_Q_vec() #Returns Q_values. If unknown, returns nan
    bus_num_init = power_network.get_bus_num_vec()  #Returns bus number as a vector (Not intitial)
    V = power_network.get_V_vec() #Returns V_values. If unknown, returns nan
    delta = power_network.get_delta_vec() #Returns delta_values. If unknown, returns nan
    PQ_vec, PQ_jacobian = power_network.get_PQ_vec() #Returns PQ_vec. This is a vector consisting of the kwnown values of P and Q. P_values first, then Q_values. (LHS of jacobi-equation) 
                                                     #Returns PQ_jacobian. Tis is av vactor of "PQ"-objects which contains information if element is P or Q, and what kind of bus the P and Q bellongs to. 
    VD_vec, VD_jacobian = power_network.get_VD_jacobian() #Returns VD_vec. This vector consists of unknown elements of V and delta. This elemens are given flat start values. 1 if V, 0 if delta. This is only valid before the first iteration. Delta first, then V. 
                                                          #Returns VD:jacobian. This vector consists of "VD"-objects which contains information if element is V or delta, and what kind of bus the V and delta bellongs to.
    bus_type_init = power_network.get_bus_type_vec() #Returns vector with bus type. O if slack, 1 if PV and 2 if PQ.       
    num_buses = len(bus_num_init)  #Calculates how many busses there are in the network. 
    
    
    delta_vd = [1] * len(VD_vec) #Initial state of delta_vd. To avoid convergens in first iteration. 
    i= 0
    while(abs(max(np.real(delta_vd))) > convergence):
        if (i==0):  #First iteration
            print("Iteration", i+1, ": \n")
            bus_type = power_network.get_bus_type_vec()
            delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V  = iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_init, delta_init, Ybus, bus_num_init, P_init, Q_init, VD_vec, power_network, bus_type_init, Q_max)
            i += 1
            
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type)
        elif (i==5):
            break
        else:
            print("Iteration", i+1, ": \n")
            delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V  = iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_updated, delta_updated, Ybus, bus_num_init, P_calc, Q_calc, VD_vec_current,power_network, bus_type, Q_max)            
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type)
            i += 1
    return P_updated, Q_updated 
      


a, b = NR(Ybus, power_network, convergence, Q_max)

print(a)
print(b)



