# In this file, the initial times follow Gamma distributions, and the sampled times during the simulation also.
import numpy as np
import matplotlib.pyplot as plt


def logdelay():
    r_p = 1 / 0.05
    k_r = 20
    k_y = 10
    k_g = 3
    cc_r = 1.028
    cc_y = 0.387
    cc_g = 1.285
    m0 = 0
    r0 = 306
    y0 = 109
    g0 = 76
    repeat = 20
    # cc_length = 10
    K = 10000
    time_step = 0.25
    Data = []
    #############################################

    for i in range(repeat):
        M = m0
        R = np.random.gamma(shape=9.223, scale=0.61, size=r0)
        Y = np.random.gamma(shape=11.777, scale=0.479, size=y0)
        G = np.random.gamma(shape=6.527, scale=0.99, size=g0)

        t = 0
        t_stat = [0]
        R_stat, Y_stat, G_stat, M_stat = [len(R)], [len(Y)], [len(G)], [M]

        while t < 125 and np.any([0 < len(R), 0 < len(Y), 0 < len(G), M < K]):  # 47.75

            tau1 = R[np.argmin(R)] if 0 < len(R) else np.Inf
            tau2 = Y[np.argmin(Y)] if 0 < len(Y) else np.Inf
            tau3 = G[np.argmin(G)] if 0 < len(G) else np.Inf
            tau4 = np.random.exponential(K / (M * r_p * (K - M - 2 * (len(R) + len(Y) + len(G))))) \
                if 0 < M * (K - M - 2 * (len(R) + len(Y) + len(G))) else np.Inf

            ####################
            min_tau = np.min(np.array([tau1, tau2, tau3, tau4]))
            t += min_tau

            if min_tau == tau1:  # Ha a legrövidebb idő R-beli sejthez tartozik
                R = np.delete(R, np.argmin(R))
                time_in_Y = np.random.gamma(k_y, cc_y)  # a sejtciklus Y-ben töltött hossza (jöhet eloszlásból)
                Y = np.append(Y, [time_in_Y])
                R -= min_tau
                Y -= min_tau
                G -= min_tau

            elif min_tau == tau2:  # Ha a legrövidebb idő Y-beli sejthez tartozik
                Y = np.delete(Y, np.argmin(Y))
                time_in_G = np.random.gamma(k_g, cc_g)  # a sejtciklus G-ben töltött hossza (jöhet eloszlásból)
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
                M -= 1
                R -= min_tau
                Y -= min_tau
                G -= min_tau
                time_in_R = np.random.gamma(k_r, cc_r)  # a sejtciklus hossza
                R = np.append(R, [time_in_R])

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

    return [avg_t, avg_R, avg_Y, avg_G, avg_M]


t, r, y, g, m = logdelay()

fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_ylabel('Nr. of cells', color="black")
ax1.plot(t, m, c="purple", linewidth=1.5, linestyle='dotted', label='Nr. of cells in M')
ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
ax1.plot(t, m + r + y + g, c='black', linewidth=2, label='Total nr. of cells')
ax1.tick_params(axis='y', labelcolor="black")
ax1.set_xlabel('Time')
plt.legend(loc='upper left')
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.savefig('changed_shape_params.png', dpi=300)
plt.show()
