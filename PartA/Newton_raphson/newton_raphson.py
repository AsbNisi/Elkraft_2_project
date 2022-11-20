import numpy as np

from Newton_raphson.NR_functions import Ybus, read_buses, iterate_NR, printing_buses, printing_Y_bus, printing_lines, Q_violated, Q_max_violation, VD_vec_Qmax, updateVD_vec, updateVD, Q_calc_violated, PQ_to_PV
from Newton_raphson.NR_network import Network

bus_vec = read_buses('PartA/Busdata.csv')

power_network = Network(bus_vec)

Ybus = Ybus('PartA/impedances.csv', len(bus_vec))

convergence = 0.00001
Q_max = [0.5,5,-1.5,5,5]

def NR(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method):
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
    bus_type_init_clean = power_network.get_bus_type_vec() #Returns vector with bus type. O if slack, 1 if PV and 2 if PQ.    
    num_buses = len(bus_num_init)  #Calculates how many busses there are in the network. 
    
    
    delta_vd = [1] * len(VD_vec) #Initial state of delta_vd. To avoid convergens in first iteration. 
    i= 0
    printing_Y_bus(Ybus)
    while(abs(max(np.real(delta_vd))) > convergence):
        if (i==0):  #First iteration
            print("Iteration", i+1, ": \n")
            bus_type = power_network.get_bus_type_vec()
            delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V  = iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_init, delta_init, Ybus, bus_num_init, P_init, Q_init, VD_vec, power_network, bus_type_init, Q_max, Q_limit, reactive_limits_method)
            i += 1
            
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type)
        else:
            print("Iteration", i+1, ": \n")
            delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V  = iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_updated, delta_updated, Ybus, bus_num_init, P_calc, Q_calc, VD_vec_current,power_network, bus_type, Q_max, Q_limit, reactive_limits_method)            
            printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type)
            i += 1

    #Her testes Q_violation etter endt iterasjon. Da må ikke Q_violated være aktivert inne i while-loopen over. 
    if(reactive_limits_method== "after"):
        if (Q_violated(Q_max, Q_updated, bus_type)):
                Q_updated, power_network = Q_max_violation(Q_updated, Q_max, bus_num_init, V, power_network)
                power_network = PQ_to_PV(bus_type_init, bus_type, power_network, V_updated)
                bus_type = power_network.get_bus_type_vec()
                VD_vec, VD_jacobian = power_network.get_VD_jacobian()
                PQ_vec, PQ_jacobian = power_network.get_PQ_vec()
                
                
                VD_vec_current = VD_vec_Qmax(VD_vec, VD_vec_current, bus_type, bus_type_init, V)
                
                
                V  = power_network.get_V_vec()
                VD_vec_current = updateVD_vec(VD_vec_current, delta_vd, bus_type_init, bus_type, delta, V)


                delta_updated, V_updated = updateVD(VD_vec_current,delta, V, bus_type_init, bus_type)
                Q_calc = Q_calc_violated(bus_type_init,bus_type, Q_updated, Q_calc)
        delta_vd = [1] * len(VD_vec) 
        while(abs(max(np.real(delta_vd))) > convergence):
                print("Iteration", i+1, ": \n")
                delta_updated, V_updated, VD_vec_current, P_calc, Q_calc, P_updated, Q_updated, bus_type, power_network, VD_jacobian, PQ_jacobian, PQ_vec, bus_type, delta_vd, V  = iterate_NR(VD_jacobian, PQ_jacobian, PQ_vec, num_buses, V, delta, V_updated, delta_updated, Ybus, bus_num_init, P_calc, Q_calc, VD_vec_current,power_network, bus_type, Q_max, Q_limit, reactive_limits_method)            
                printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type)
                i += 1

    Power_network = PQ_to_PV(bus_type_init_clean, bus_type, power_network, V_updated) #Sets the transfrormed PV_bus back to a PV_bus.
    printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_init_clean)
    printing_lines(bus_type_init_clean, "PartA/impedances.csv", V_updated, Ybus, delta_updated)
    return P_updated, Q_updated 

    
      



