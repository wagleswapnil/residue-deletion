import numpy as np
rAmax = 0.2
rAmin = 0.2

rBmax = 0.1334
rBmin = 0.1334

fcA = 41003.20000
fcB = 41003.20000

bonded_l ="0.00 0.10 0.20 0.30 0.40 0.50 0.60 0.65 0.70 0.75 0.80 0.85 0.90 0.95 0.99 1.00"
arr = [float(x) for x in bonded_l.split()]
bonded_lambdas = np.array(arr)
restraint_lambdas = np.array([0.00, 0.00, 0.10, 0.20, 0.30, 0.30, 0.40, 0.50, 0.60, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])


for idx, lambda_element in enumerate(bonded_lambdas):
    rmax = rAmax * (1 - lambda_element) + rBmax * lambda_element
    rmin = rAmin * (1 - lambda_element) + rBmax * lambda_element
    fc = fcA * (1 - lambda_element) + fcB * lambda_element
    print(f"idx={idx:3d}, lambda={lambda_element:5.2f}, rmax={round(rmax, 4):7.4f}, rmin={round(rmin, 4):7.4f}, fc={round(fc, 2):12.7f}")
