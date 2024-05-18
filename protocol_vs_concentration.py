import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from protocol_for_experiment import logdelay
import time

start = time.time()
parameters_1 = [list(np.linspace(0,120,2)), list(np.linspace(0,120,4)), list(np.linspace(0,120,6)), list(np.linspace(0,120,8)), list(np.linspace(0,120,10)), list(np.linspace(0,120,12)), list(np.linspace(0,120,14)), list(np.linspace(0,120,16)), list(np.linspace(0,120,18)), list(np.linspace(0,120,20))]  # adagolások ideje

hm_data = np.zeros(len(parameters_1))

for i in range(len(parameters_1)):  # sorokban protocol
    t, r, y, g, m, d = logdelay(protocol=parameters_1[i])
    total_unhealty = r + y + g + m
    hm_data[i] = total_unhealty[-1] / total_unhealty[0]
    print(i)

end = time.time()
print('One loop:', end - start)

treatments = [str(2), str(4), str(6), str(8), str(10), str(12), str(14), str(16), str(18), str(20)]
values = hm_data

fig = plt.figure(figsize=(10, 5))

# creating the bar plot
plt.bar(treatments, values, color='maroon')

plt.xlabel("Number of treatments")
plt.ylabel("Relative cell count")
plt.title("Relative cell count across different dosage strategies")
plt.savefig('dose_strat_3', dpi=300)
plt.show()
