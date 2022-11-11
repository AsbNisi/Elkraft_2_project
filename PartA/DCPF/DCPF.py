import numpy as np
from Newton_raphson.NR_functions import read_buses, P_Calc, Ybus
from DCLF.DCLF_functions import Ybus_dclf
from Newton_raphson.NR_network import Network
from FDCLF.FDCLF_functions import Ybus_fdclf
import timeit
power_network = Network(read_buses('PartA/Busdata.csv'))
Ybus = Ybus('PartA/impedances.csv', 5)

# Attaining the DCPF Ybus for a given power network
def Ybus_DCPF(power_network):
    BusNum =len(power_network.buses)
    
    # B' matrix from FDLF can be multiplied by -1 to attain DCPF Ybus
    b_dash = Ybus_fdclf('PartA/impedances.csv', BusNum, 'PartA/Busdata.csv', power_network, Ybus)[0]
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
    Ybus_DCLF = Ybus_dclf('PartA/impedances.csv', BusNum)
    P_injections = np.real(P_Calc(V_o, Ybus_DCLF, range(BusNum), delta_vec, P_o))
    
    return P_injections, delta_vec


start_time = timeit.default_timer()

power_network = Network(read_buses('PartA/Busdata.csv'))
P_injections, phase_angles = DCPF_calc(power_network)

runtime = timeit.default_timer() - start_time
print(f'Runtime: {runtime}')  
