import os
import sys


w_fnet = float(sys.argv[1])
w_pnet = float(sys.argv[2])

S_fnet = int(sys.argv[3])
S_pnet = int(sys.argv[4])

T2_fnet = float(sys.argv[5])
T2_pnet = float(sys.argv[6])

lambd = float(sys.argv[7])

rand_seed = 1

Z_flag = int(sys.argv[8])
Z = int(sys.argv[9])

print("Starting simulation!")
os.system(
    f"./backend/model_potts_hybrid/runsim.exe {w_pnet} {w_fnet} {S_pnet} {S_fnet} {T2_pnet} {T2_fnet} {lambd} {lambd} {rand_seed} {Z_flag} {Z}"
)
print("Yaay!")
