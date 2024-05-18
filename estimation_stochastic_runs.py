import numpy as np
import matplotlib.pyplot as plt
from data import *


def logdelay(params):
    r_p, k_r, cc_r, k_y, cc_y, k_g, cc_g, k_r0, cc_r0, k_y0, cc_y0, k_g0, cc_g0 = params
    m0 = 0
    r0 = 106
    y0 = 209
    g0 = 76
    repeat = 200
    K = 10000
    time_step = 0.1
    Data = []
    Data_raw = []
    Distributions = []
    #############################################

    for i in range(repeat):
        M = m0

        R = np.random.gamma(k_r0, cc_r0, r0)  # uniform(0, 2*10*1.028, r0)
        Y = np.random.gamma(k_y0, cc_y0, y0)  # uniform(0, 2*10*0.387, y0)
        G = np.random.gamma(k_g0, cc_g0, g0)  # np.random.gamma(k_g0, cc_g0, g0) #

        t = 0
        t_stat = [0]
        D_raw = [[t, len(R), len(Y), len(G), M]]
        R_stat, Y_stat, G_stat, M_stat = [len(R)], [len(Y)], [len(G)], [M]
        # 9.75
        while t < 50 + time_step and np.any([0 < len(R), 0 < len(Y), 0 < len(G), M < K]):

            tau1 = R[np.argmin(R)] if 0 < len(R) else np.Inf
            tau2 = Y[np.argmin(Y)] if 0 < len(Y) else np.Inf
            tau3 = G[np.argmin(G)] if 0 < len(G) else np.Inf
            tau4 = np.random.exponential(K / (M * r_p * (K - M - 2 * (len(R) + len(Y) + len(G))))) \
                if 0 < M * (K - M - 2 * (len(R) + len(Y) + len(G))) else np.Inf

            ####################
            min_tau = np.min([tau1, tau2, tau3, tau4])
            t += min_tau

            if min_tau == tau1:  # Ha a legrövidebb idő R-beli sejthez tartozik
                R = np.delete(R, np.argmin(R))
                time_in_Y = np.random.gamma(k_y,
                                            cc_y) + min_tau  # a sejtciklus Y-ben töltött hossza (jöhet eloszlásból)
                Y = np.append(Y, [time_in_Y])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau2:  # Ha a legrövidebb idő Y-beli sejthez tartozik
                Y = np.delete(Y, np.argmin(Y))
                time_in_G = np.random.gamma(k_g,
                                            cc_g) + min_tau  # a sejtciklus G-ben töltött hossza (jöhet eloszlásból)
                G = np.append(G, [time_in_G])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau3:  # Ha a legrövidebb idő G-beli sejthez tartozik
                G = np.delete(G, np.argmin(G))
                M += 2
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau4:  # Enter cell cycle
                time_in_R = np.random.gamma(k_r, cc_r) + min_tau  # a sejtciklus hossza
                R = np.append(R, [time_in_R])
                M -= 1
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            D_raw.append([t, len(R), len(Y), len(G), M])
            ######################## Statistics ################################

            disc_step = int((t - t_stat[-1]) / time_step)

            if disc_step >= 0:
                for j in range(disc_step):
                    t_stat.append(t_stat[-1] + time_step)
                    R_stat.append(len(R))
                    Y_stat.append(len(Y))
                    G_stat.append(len(G))
                    M_stat.append(M)

        ####################################################################
        Data.append([t_stat, R_stat, Y_stat, G_stat, M_stat])
        D_raw = np.array(D_raw).T
        Data_raw.append(D_raw)
        Distributions.append([R, Y, G])

    m = len(Data[0][0])
    size = len(Data)
    for d in Data:
        if len(d[0]) < m:
            m = len(d[0])

    avg_t, avg_R, avg_Y, avg_G, avg_M = np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m)

    for d in Data:
        avg_t += np.array(d[0])[:m] / size
        avg_R += np.array(d[1])[:m] / size
        avg_Y += np.array(d[2])[:m] / size
        avg_G += np.array(d[3])[:m] / size
        avg_M += np.array(d[4])[:m] / size

    return [avg_t, avg_R, avg_Y, avg_G, avg_M, Data_raw, Distributions]


pars = [5.101089691028263e+01, 4.134421316370464e+00, 1.460202315290193e+00, 1.253527743481907e+01,
        7.430237421109924e-01, 8.418714315835516e+00, 6.563746706018511e-01, 8.947126588628151e+00,
        7.010408140809474e-01, 1.163077742270154e+01, 4.571958301088177e-01, 5.861721344080051e+00,
        9.650462027760387e-01]

t, r, y, g, m, data, distributions = logdelay(pars)

if True:
    for i in range(len(data)):
        t_stat, R_stat, Y_stat, G_stat, M_stat = data[i]
        # plt.plot(t_stat, np.array(M_stat), c="black", linewidth=1, alpha = 0.1)
        plt.plot(t_stat, np.array(M_stat) + np.array(R_stat), c="firebrick", linewidth=1, alpha=0.1, label='Simulated R(t)')
        plt.plot(t_stat, Y_stat, c="darkorange", linewidth=1, alpha=0.1, label='Simulated Y(t)')
        plt.plot(t_stat, G_stat, c="darkgreen", linewidth=1, alpha=0.1, label='Simulated G(t)')

plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[0], linewidth=2, alpha=0.9, c="mistyrose", label="Measured no. of cells in R")
plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[1], linewidth=2, alpha=0.9, c="yellow", label="Measured no. of cells in Y")
plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[2], linewidth=2, alpha=0.9, c="mediumseagreen", label="Measured no. of cells in G")
plt.title('The number of cells in the different cell cycle phases across time')
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.legend()
plt.savefig('estimated_stoch_runs2', dpi=300)
plt.show()
