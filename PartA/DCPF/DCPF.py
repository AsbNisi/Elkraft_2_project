import numpy as np
from Newton_raphson.NR_functions import read_buses, P_Calc, Ybus, printing_buses
from DCLF.DCLF_functions import Ybus_dclf
from Newton_raphson.NR_network import Network
from FDCLF.FDCLF_functions import Ybus_fdclf
import pandas as pd
import cmath
bus_vec = read_buses('Part1/Busdata.csv')
power_network = Network(bus_vec)
Ybus = Ybus('Part1/impedances.csv', 5)

# Attaining the DCPF Ybus for a given power network
def Ybus_DCPF(power_network):
    BusNum =len(power_network.buses)
    
    bus_type_vec = power_network.get_bus_type_vec()
    # B' matrix from FDLF can be multiplied by -1 to attain DCPF Ybus
    b_dash = Ybus_fdclf(bus_type_vec, BusNum, Ybus)[0]
    Y_dcpf = np.real(b_dash.copy() * -1)
    return Y_dcpf


# Running DCPF study on a given power network
def DCPF_calc(power_network):
    BusNum =len(power_network.buses)  
    
    # Ybus for DCPF is attained and inversed
    Y_dcpf = Ybus_DCPF(power_network)
    Y_dcpf_inv = np.linalg.inv(Y_dcpf)
    
    # Column vector of powers is attained and the slack bus removed
    P_vec = power_network.get_P_vec()
    slack_bus_index = np.where(np.isnan(P_vec))[0]
    del P_vec[slack_bus_index[0]]
    P_vec = np.array([P_vec]).T
    
    # Power angles are attained by DCPF, then the slack bus is reintroduced
    delta_vec_reduced = np.matmul(Y_dcpf_inv, P_vec)
    delta_vec = np.insert(delta_vec_reduced, slack_bus_index[0], 0)
    
    # For DCPF, Q values are omitted, and all voltages are 1pu
    P_o = np.array([0 for _ in range(BusNum)])
    V_o = np.array([1 for _ in range(BusNum)])
    
    # Ybus for DCLF is used to find power flows as it only includes line susceptances
    Ybus_DCLF = Ybus_dclf('impedances.csv', BusNum)
    P_injections = np.real(P_Calc(V_o, Ybus_DCLF, range(BusNum), delta_vec, P_o))
    
    bus_num_init = power_network.get_bus_num_vec()
    bus_type_init_clean = power_network.get_bus_type_vec()
    P_updated = P_injections
    Q_updated = [0 for _ in range(len(P_updated))]
    V_updated = [1 for _ in range(len(P_updated))]
    delta_updated = delta_vec
    printing_buses(V_updated, delta_updated, P_updated, Q_updated, bus_num_init, bus_type_init_clean)
    printing_lines(bus_vec, "impedances.csv", V_updated, Ybus, delta_updated)
    return P_injections, delta_vec



#Printing load flow in lines
def printing_lines(bus_vec, file, V_updated, Ybus, delta_updated):
    
    #Calculate complex power flow in lines
    df_lines = pd.read_csv(file, sep=";")
    V_updated = np.real(V_updated)
    delta_updated = np.real(delta_updated)
    
    V_complex = np.zeros(len(bus_vec), dtype=complex)
    for i in range(len(bus_vec)):
        V_complex[i] = cmath.rect(V_updated[i], delta_updated[i])
    
    S_ik = np.zeros((len(bus_vec), len(bus_vec)), dtype=complex)
    for i in range(len(bus_vec)):
        shunt = df_lines["Full_line_B"][i]/2*(1j)
        
        for k in range(len(bus_vec)):
            
            S_ik[i][k] = V_complex[i]*np.conjugate(-Ybus[i][k]*(V_complex[i]-V_complex[k]) + shunt*V_complex[i])#*S_base
            
    print("Updated line info:", '\n')
    
    d = {}
    for i in range (len(bus_vec)):
        from_line = df_lines["From_line"][i]
        to_line = df_lines["To_line"][i]

        line = str(from_line) + " - " + str(to_line)
        
        d[line] = S_ik[from_line-1,to_line-1], np.real(S_ik[from_line-1,to_line-1]), np.imag(S_ik[from_line-1,to_line-1]), abs(np.real(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1])), abs(np.imag(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1])) 
    print ("{:<7} {:<12} {:<12} {:<12} {:<12}".format('Line','Active Power Flow [MW] ', 'Reactive Power Flow [MVar] ', "Active Power Loss [MW] ", "Reactive Power Loss [MVar] "))    
    for k, v in d.items():
        apparent, active, reactive, ploss, qloss = v
        print("{:<7} {:<23} {:<27} {:<23} {:<19}".format(k, round(active,4), 0, round(ploss,4), 0))
    print('\n')
    
    
    e = {}
    for i in range (len(bus_vec)):
        to_line = df_lines["From_line"][i]
        from_line = df_lines["To_line"][i]

        line = str(from_line) + " - " + str(to_line)
        
        e[line] = S_ik[from_line-1,to_line-1], np.real(S_ik[from_line-1,to_line-1]), np.imag(S_ik[from_line-1,to_line-1])#, np.real(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1]), np.imag(S_ik[from_line-1,to_line-1]+S_ik[to_line-1,from_line-1]) 
    for k, v in e.items():
        apparent, active, reactive = v
        print("{:<7} {:<23} {:<27}".format(k, round(active,4), 0))#, round(ploss,4), round(qloss,4)))
    print('\n')
    
    return 