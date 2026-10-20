import os, sys
import decimal
in_file = 'production.log'
in_replica = sys.argv[1].strip()
current_idx = in_replica
out_file_replica = open("out_file_replica.csv", 'w')
with open(in_file, 'r') as f:
    for line in f:
        if "Replica exchange at step" in line:
            t = line.split()[-1]
            line = f.readline()
            if "dE_term =" in line:
                line = f.readline()

            info = line.split()
            try:
                idx = info.index(current_idx)
            except:
                print("something went wrong")
                quit()
            if current_idx != "14" and info[idx + 1] == "x":
                current_idx = info[idx + 2]
            elif info[idx - 1] == "x":
                current_idx = info[idx - 2]
            else:
                pass
            out_file_replica.write(t + "  " + current_idx + "\n")


out_file_replica.close()


