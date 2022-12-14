from Newton_raphson.newton_raphson import NR, power_network, Ybus
from DCLF.decoupled_load_flow import DCLF, power_network
from FDCLF.fast_decoupled_load_flow import FDCLF
from DCPF.DCPF import DCPF_calc
from Task4.Task4_NR import NR_trans, Ybus_trans, power_network_trans
from Task4.Task4_FDCLF import FDCLF_trans
import time


# Setting convergence criteria
convergence = 0.000000001

#Vector with reactive power limits for the base case
Q_max = [0.5,5,-1.5,5,5]
#Vector with reactive power limits for the case with transformers
Q_max_trans = [0.5,5,-1.5,5,5, 5, 5]


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
          If DC power flow write 'DCPF'\n
          If Newton-Raphson with tansformers write 'NR_trans',\n
          If Fast Decoupled method with transformers write 'FDCLF_trans'\n """)
    
    reactive_limits = input("Do you want to run the Load Flow Analysis with reactive power limits: 'y'/'n' \n")
    if(reactive_limits == 'y'):
        reactive_limits_method = input("Do you want to check Q_violation before or after ended iteration? 'before'/'after' \n")
    else:
        reactive_limits_method = None
    
    print("Chosen method: ", method, "\n\n")
    
    if (method == "NR" and reactive_limits == "n"):
        Q_limit = False
        start = time.process_time()
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == "NR" and reactive_limits == "y"):
        Q_limit = True
        start = time.process_time()
        P_updated, Q_updated = NR(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'DCLF' and reactive_limits == "n"):
        Q_limit = False
        start = time.process_time()
        P_updated, Q_updated  = DCLF(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'DCLF' and reactive_limits == "y"):
        Q_limit = True
        start = time.process_time()
        P_updated, Q_updated  = DCLF(Ybus, power_network, convergence, Q_max, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'FDCLF_1' and reactive_limits == "n"):
        method = 1
        Q_limit = False
        start = time.process_time()
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'FDCLF_1' and reactive_limits == "y"):
        method = 1
        Q_limit = True
        start = time.process_time()
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)   
        print(time.process_time() - start)
    if (method == 'FDCLF_2' and reactive_limits == "n"):
        method = 2
        Q_limit = False
        start = time.process_time()
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'FDCLF_2' and reactive_limits == "y"):
        method = 2
        Q_limit = True
        start = time.process_time()
        P_updated, Q_updated = FDCLF(Ybus, power_network, convergence, Q_max, method, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'DCPF'):
        start = time.process_time()
        P_injections, delta_vec = DCPF_calc(power_network)
        print(time.process_time() - start)
    if (method == "NR_trans" and reactive_limits == "n"):
        Q_limit = False
        start = time.process_time()
        P_updated, Q_updated = NR_trans(Ybus_trans, power_network_trans, convergence, Q_max_trans, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == "NR_trans" and reactive_limits == "y"):
        Q_limit = True
        start = time.process_time()
        P_updated, Q_updated = NR_trans(Ybus_trans, power_network_trans, convergence, Q_max_trans, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'FDCLF_trans' and reactive_limits == "n"):
        Q_limit = False
        method = 1
        start = time.process_time()
        P_updated, Q_updated = FDCLF_trans(Ybus_trans, power_network_trans, convergence, Q_max_trans, method, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    if (method == 'FDCLF_trans' and reactive_limits == "y"):
        Q_limit = True
        method = 1
        start = time.process_time()
        P_updated, Q_updated = FDCLF_trans(Ybus_trans, power_network_trans, convergence, Q_max_trans, method, Q_limit, reactive_limits_method)
        print(time.process_time() - start)
    

main()