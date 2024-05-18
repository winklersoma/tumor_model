import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from logdelay import logdelay

parameters_1 = [3, 5, 7, 9, 11, 13, 15]  # adagolások száma / sűrűsége
parameters_2 = [0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.4]  # D_0 kezdeti mennyiség

hm_data = np.zeros((len(parameters_1), len(parameters_2)))

for i in range(len(parameters_1)):  # sorokban treatment frequency
    for j in range(len(parameters_2)):
        t, r, y, g, m, d = logdelay(treatment_freq=parameters_1[i], drug_0=parameters_2[j])
        total_unhealty = r + y + g + m
        hm_data[i][j] = total_unhealty[-1] / total_unhealty[0]

df = pd.DataFrame(hm_data, index=parameters_1, columns=parameters_2)
hm = sns.heatmap(data=df, annot=True)
plt.xlabel('Dosage size (in mg)')
plt.ylabel('Treatment frequency (hours between two treatments)')
plt.title("Relative count of unhealty cells")
plt.show()


fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_ylabel('Nr. of cells', color="black")
ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
ax1.plot(t, r+y+g, c='black', linewidth=2, label='Total nr. of cells')
ax1.tick_params(axis='y', labelcolor="black")
ax1.set_xlabel('Time')

ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-ax1
ax2.set_ylabel('Drug level', color="blue")
ax2.step(t, d, "blue", where='pre', alpha=0.4)
ax2.tick_params(axis='y', labelcolor="blue")
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.show()
