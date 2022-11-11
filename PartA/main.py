from Newton_raphson.newton_raphson import NR, power_network, Ybus, convergence, Q_max
from DCLF.decoupled_load_flow import DCLF, power_network, Ybus_dclf, convergence, Q_max
from FDCLF.fast_decoupled_load_flow import FDCLF
from DCPF import DCPF
Ybus_dclf = Ybus_dclf('PartA/impedances.csv', 5)

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
          If Fast Decoupled write 'FDCLF',
          If DC power flow write 'DCPF'\n""")
    
    reactive_limits = input("Do you want to run the Load Flow Analysis with reactive power limits: y/n \n")
    
    print("Chosen method: ", method, "\n\n")
    
    if (method == "NR" and reactive_limits == "n"):
        Q_limit = False
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit)
    if (method == "NR" and reactive_limits == "y"):
        Q_limit = True
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit)
    if (method == 'DCLF' and reactive_limits == "n"):
        Q_limit = False
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit)
    if (method == 'DCLF' and reactive_limits == "y"):
        Q_limit = True
        P_updated, Q_updated  = DCLF(Ybus_dclf, power_network, convergence, Q_max, Q_limit)
    if (method == 'FDCLF'):
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max)
    if (method == 'DCPF'):
        P_injections, delta_vec = DCPF(power_network)
    

main()
