from Newton_raphson.newton_raphson import NR, power_network, Ybus
from DCLF.decoupled_load_flow import DCLF, power_network, Ybus_dclf
from FDCLF.fast_decoupled_load_flow import FDCLF
from DCPF.DCPF import DCPF_calc
from Task4.Task4_NR import NR_trans, Ybus_trans, power_network_trans
from Task4.Task4_FDCLF import FDCLF_trans

Ybus_dclf = Ybus_dclf('PartA/impedances.csv', 5)
convergence = 0.000000001
Q_max = [0.5,5,-1.5,5,5]

#-------------------------------------------------------------------------
# Instructions
# Choose the method you want to use to solve the load flow problem
# A methos is chosen by typing "'method-name'" into the main function
# The following methods can be chosen
# a) "NR" - Newton Raphson
# b) "NR_Q_max" - Newton Raphson with Q_max limits included
# c) "DCLF" - Decoupled load flow
# d) "DCLF_Q_max" - Decoupled load flow with Q_max limits included
# e) FDCLF - Fast decoupled load flow
# f) DCPF - DC power flow
#-------------------------------------------------------------------------


def main():
    method = input("""Which method do you want to run? \n
          If Newton-Raphson write 'NR',\n
          If Decoupled method write 'DCLF',\n
          If Fast Decoupled method 1 (update the partial initial estimates half-way through the algorithm) write 'FDCLF_1', \n
          If Fast Decoupled method 2 (update the initial estimates only at the end of the first iteration) write 'FDCLF_2', \n
          If DC power flow write 'DCPF'\n""")
    
    reactive_limits = input("Do you want to run the Load Flow Analysis with reactive power limits: 'y'/'n' \n")
    if(reactive_limits == 'y'):
        reactive_limits_method = input("Do you want to check Q_violation before or after ended iteration? 'before'/'after' \n")
    else:
        reactive_limits_method = None
    
    print("Chosen method: ", method, "\n\n")
    
    if (method == "NR" and reactive_limits == "n"):
        Q_limit = False
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == "NR" and reactive_limits == "y"):
        Q_limit = True
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == 'DCLF' and reactive_limits == "n"):
        Q_limit = False
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == 'DCLF' and reactive_limits == "y"):
        Q_limit = True
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_1' and reactive_limits == "n"):
        method = 1
        Q_limit = False
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_1' and reactive_limits == "y"):
        method = 1
        Q_limit = True
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_2' and reactive_limits == "n"):
        method = 2
        Q_limit = False
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_2' and reactive_limits == "y"):
        method = 2
        Q_limit = True
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
    if (method == 'DCPF'):
        P_injections, delta_vec = DCPF_calc(power_network)
    if (method == "NR_trans" and reactive_limits == "n"):
        Q_limit = False
        P_updated, Q_updated = NR_trans(Ybus_trans, power_network_trans, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == "NR_trans" and reactive_limits == "y"):
        Q_limit = True
        P_updated, Q_updated = NR_trans(Ybus_trans, power_network_trans, convergence, Q_max, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_trans' and reactive_limits == "n"):
        Q_limit = False
        method = 1
        P_updated, Q_updated = FDCLF_trans(Ybus_trans, power_network_trans, convergence, Q_max, method, Q_limit, reactive_limits_method)
    if (method == 'FDCLF_trans' and reactive_limits == "y"):
        Q_limit = True
        method = 1
        P_updated, Q_updated = FDCLF_trans(Ybus_trans, power_network_trans, convergence, Q_max, method, Q_limit, reactive_limits_method)

    

main()