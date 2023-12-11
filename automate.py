# print("hello world 2")
import os
import sys


S = int(sys.argv[1])
w = float(sys.argv[2])
T2 = int(sys.argv[3])

print("Starting simulation!")
os.system(f"backend/model_potts_homo/runsim {S} {w} {T2}")
print("Yaay!")
