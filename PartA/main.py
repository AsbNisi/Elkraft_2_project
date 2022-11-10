from Newton_raphson.newton_raphson import NR, power_network, Ybus, convergence, Q_max
from DCLF.decoupled_load_flow import DCLF, power_network, Ybus_dclf, convergence, Q_max
from FDCLF.fast_decoupled_load_flow import FDCLF
from DCPF import DCPF
from Task4.Task4_NR import NR_trans, power_network, Ybus, convergence, Q_max
from Task4.Task4_FDCLF import FDCLF_trans, power_network, Ybus, convergence, Q_max


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
          If Newton-Raphson write 'NR', with reactive limits write 'NR_Q_max',\n
          If Decoupled method write 'DCLF', with reactive limits write 'DCLF_Q_max',\n
          If Fast Decoupled write 'FDCLF',
          If DC power flow write 'DCPF'\n""")
    
    print("Chosen method: ", method, "\n\n")
    
    if (method == "NR"):
        Q_limit = False
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit)
    if (method == "NR_Q_max"):
        Q_limit = True
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit)
    if (method == 'DCLF'):
        Q_limit = False
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit)
    if (method == 'DCLF_Q_max'):
        Q_limit = True
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit)
    if (method == 'FDCLF'):
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max)
    if (method == 'DCPF'):
        P_injections, delta_vec = DCPF(power_network)
    if (method == 'NR_trans'):
        Q_limit = False
        P_updated, Q_updated = NR_trans(Ybus, power_network, convergence, Q_max, Q_limit)
    if (method == 'FDCLF_trans'):
        Q_limit = False
        P_updated, Q_updated = FDCLF_trans(Ybus, power_network, convergence, Q_max, Q_limit)
    

main()
