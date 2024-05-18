from data import *
from functions import cc_3_logdelay

plt.scatter(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[0], linewidth=0.5, alpha=0.4, c="red")
plt.scatter(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[1], linewidth=0.5, alpha=0.4, c="orange")
plt.scatter(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[2], linewidth=0.5, alpha=0.4, c="green")
plt.title('The number of cells in the different cell cycle phases across time')
plt.xlabel('Time')
plt.ylabel('Number of cells')

t, r, y, g, m = cc_3_logdelay(
    params=[5.101089691028263e+01, 4.134421316370464e+00, 1.460202315290193e+00, 1.253527743481907e+01,
            7.430237421109924e-01, 8.418714315835516e+00, 6.563746706018511e-01, 8.947126588628151e+00,
            7.010408140809474e-01, 1.163077742270154e+01, 4.571958301088177e-01, 5.861721344080051e+00,
            9.650462027760387e-01])

plt.plot(t, m + r, c="red", linewidth=1.5, label='Nr. of cells in R')
plt.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
plt.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
# plt.savefig('estimated_vs_real.png', dpi=300)
plt.show()
