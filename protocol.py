import numpy as np
import matplotlib.pyplot as plt
import random
from data import *
import time

start = time.time()


def logdelay(drug_0=2.3869, k_c=(0, 0, 0.9), delta=1.8328, protocol=[0, 24, 48, 72, 96], half_life=2):
    r_p = 5.101089691028263e+01
    k_r = 4.134421316370464e+00
    k_y = 1.253527743481907e+01
    k_g = 8.418714315835516e+00
    cc_r = 1.460202315290193e+00
    cc_y = 7.430237421109924e-01
    cc_g = 6.563746706018511e-01
    # print("Az R osztályban töltött idő várható értéke: " + str(k_r * cc_r))
    # print("Az Y osztályban töltött idő várható értéke: " + str(k_y * cc_y))
    # print("A G osztályban töltött idő várható értéke: " + str(k_g * cc_g))
    # print("A sejtciklusban töltött idő várható értéke: " + str(k_r * cc_r + k_y * cc_y + k_g * cc_g))
    drug_decay = np.log(2) / half_life  # 30.8  # derived from Joerger et al.
    delta = delta  # drug elimination rate from de Pillis et al.
    drug_0 = drug_0
    K_c_r = k_c[0]  # Selectivity
    K_c_y = k_c[1]  # Selectivity
    K_c_g = k_c[2]  # Selectivity
    m0 = 0
    r0 = 3000
    y0 = 3000
    g0 = 3000
    repeat = 2
    K = 100000  # meg kellett emelni, mert elérte. (K - 2 * M * (len(R) + len(Y) + len(G))) < 0!!
    end_time = protocol[-1] + 24 if protocol[-1] != 0 else 120
    time_step = 0.05
    drug_step = 0.05
    Data = []

    #############################################

    for i in range(repeat):

        if protocol[0] == 0:
            D = drug_0
            drug_time = drug_step  # utolsó kezelés óta eltelt idő
            treat_count = 1
        else:
            D = 0
            drug_time = np.inf  # utolsó kezelés óta eltelt idő
            treat_count = 0

        if len(protocol) > treat_count:
            next_treatment = protocol[treat_count]  # mikor esedékes a következő kezelés
        else:
            next_treatment = np.inf
        """    
        D = drug_0 if protocol[0] == 0 else 0
        drug_time = drug_step if protocol[0] == 0 else np.inf  # az utolsó kezelés óta eltelt idő

        if protocol[0] == 0:
            protocol.pop(0)
        treat_count = 0  # számolja a kezeléseket
        if len(protocol) > 0:
            next_treatment = protocol[treat_count] # if protocol != [0] else np.inf  # mikor esedékes a következő kezelés
        else:
            next_treatment = np.infa
        """
        M = m0
        R = np.random.gamma(shape=8.947126588628151e+00, scale=7.010408140809474e-01, size=r0)
        Y = np.random.gamma(shape=1.163077742270154e+01, scale=4.571958301088177e-01, size=y0)
        G = np.random.gamma(shape=5.861721344080051e+00, scale=9.650462027760387e-01, size=g0)

        t = 0
        t_stat, t_drug = [0], [0]

        R_stat, Y_stat, G_stat, M_stat, D_stat = [len(R)], [len(Y)], [len(G)], [M], [D]

        while t < end_time and np.any([0 < len(R), 0 < len(Y), 0 < len(G)]) and (
                (M + 2 * (len(R) + len(Y) + len(G))) < K):
            # Gillespie
            a_r = r_p * M * (K - (M + 2 * (len(R) + len(Y) + len(G)))) / K
            a_d_r = K_c_r * len(R) * (1 - np.exp(-delta * D))
            a_d_y = K_c_y * len(Y) * (1 - np.exp(-delta * D))
            a_d_g = K_c_g * len(G) * (1 - np.exp(-delta * D))
            a = [a_r, a_d_r, a_d_y, a_d_g]
            a0 = sum(a)

            Tau = np.random.exponential(scale=(1 / a0)) if 0 < a0 else np.inf
            theta1 = R[np.argmin(R)] if 0 < len(R) else np.Inf
            theta2 = Y[np.argmin(Y)] if 0 < len(Y) else np.Inf
            theta3 = G[np.argmin(G)] if 0 < len(G) else np.Inf

            time_list = [theta1, theta2, theta3, drug_time, next_treatment]
            if Tau < min(time_list):
                min_tau = Tau

                mu = random.choices([0, 1, 2, 3], weights=a, k=1)

                if mu[0] == 0:  # a_r nyert -> újabb sejt kerül a sejtciklusba
                    M -= 1
                    time_in_R = np.random.gamma(k_r, cc_r) + min_tau  # a sejtciklus R-ben töltött véletlen hossza
                    R = np.append(R, [time_in_R])
                elif mu[0] == 1:  # a_d_r nyert
                    choice = np.random.choice(len(R), size=1)
                    R = np.delete(R, choice)
                elif mu[0] == 2:  # a_d_y nyert
                    choice = np.random.choice(len(Y), size=1)
                    Y = np.delete(Y, choice)
                else:  # a_d_g nyert
                    choice = np.random.choice(len(G), size=1)
                    G = np.delete(G, choice)

            elif min(time_list) == theta1:  # Ha a legrövidebb idő R-beli sejthez tartozik
                # print("R")
                # print("___")
                min_tau = theta1
                R = np.delete(R, np.argmin(R))
                time_in_Y = np.random.gamma(k_y,
                                            cc_y) + min_tau  # a sejtciklus Y-ben töltött hossza (jöhet eloszlásból)
                Y = np.append(Y, [time_in_Y])

            elif min(time_list) == theta2:  # Ha a legrövidebb idő Y-beli sejthez tartozik
                # print("Y")
                # print("___")
                min_tau = theta2
                Y = np.delete(Y, np.argmin(Y))
                time_in_G = np.random.gamma(k_g,
                                            cc_g) + min_tau  # a sejtciklus G-ben töltött hossza (jöhet eloszlásból)
                G = np.append(G, [time_in_G])

            elif min(time_list) == theta3:  # Ha a legrövidebb idő G-beli sejthez tartozik
                # print("G")
                # print("___")
                min_tau = theta3
                G = np.delete(G, np.argmin(G))
                M += 2

            elif min(time_list) == drug_time:  # itt van a bomlás diszkretizált időben
                # print(t + drug_time)
                # print("Drug decay")
                # print("drug_time: " + str(drug_time))
                # print("next_treatment: " + str(next_treatment))
                # print(10*"__")
                min_tau = drug_time
                drug_time = drug_step + min_tau  # ez a lépés azért van, mert lejjebb MINDENBŐL kivonjuk a min_tau-t
                D = D * np.exp(- drug_decay * drug_step)

            elif min(time_list) == next_treatment:
                # print(t + next_treatment)
                treat_count += 1
                # print("Next treatment")
                # print("drug_time: " + str(drug_time))
                # print("next_treatment: " + str(next_treatment))
                # print("treat_count: " + str(treat_count))
                # print(10*"__")
                min_tau = next_treatment
                if drug_time == np.inf:
                    drug_time = drug_step

                if treat_count < len(protocol):
                    next_treatment = protocol[treat_count] - protocol[treat_count - 1] + min_tau
                    D += drug_0
                elif drug_0 != 0:
                    # next_treatment = protocol[-1] - protocol[-2] + min_tau
                    next_treatment = np.inf
                    D += drug_0

            R -= min_tau
            Y -= min_tau
            G -= min_tau
            drug_time -= min_tau
            next_treatment -= min_tau

            t += min_tau

            ######################## Statistics ################################

            disc_step = int((t - t_stat[-1]) / time_step)

            if disc_step >= 0:
                for j in range(disc_step):
                    t_stat.append(t_stat[-1] + time_step)
                    R_stat.append(len(R))
                    Y_stat.append(len(Y))
                    G_stat.append(len(G))
                    M_stat.append(M)
                    D_stat.append(D)

        ####################################################################
        Data.append([t_stat, R_stat, Y_stat, G_stat, M_stat, D_stat])

    m = len(Data[0][0])
    size = len(Data)
    for d in Data:
        if len(d[0]) < m:
            m = len(d[0])

    avg_t, avg_R, avg_Y, avg_G, avg_M, avg_D = (
        np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m), np.zeros(m))

    for d in Data:
        avg_t += np.array(d[0])[:m] / size
        avg_R += np.array(d[1])[:m] / size
        avg_Y += np.array(d[2])[:m] / size
        avg_G += np.array(d[3])[:m] / size
        avg_M += np.array(d[4])[:m] / size
        avg_D += np.array(d[5])[:m] / size
    return [avg_t, avg_R, avg_Y, avg_G, avg_M, avg_D]


t, r, y, g, m, d = logdelay(protocol=[0, 12, 24])

end = time.time()
print('One loop:', end - start)

fig, ax1 = plt.subplots(figsize=(8, 6))

ax1.set_ylabel('Nr. of cells', color="black")
ax1.plot(t, r, c="red", linewidth=1.5, label='Nr. of cells in R')
ax1.plot(t, y, c="orange", linewidth=1.5, label='Nr. of cells in Y')
ax1.plot(t, g, c="green", linewidth=1.5, label='Nr. of cells in G')
# ax1.plot(t, r + y + g, c='black', linewidth=2, label='Total nr. of cells')
ax1.tick_params(axis='y', labelcolor="black")
ax1.set_xlabel('Time')

if any(d) != 0 or 0.0:
    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-ax1
    ax2.set_ylabel('Drug concentration, $\\frac{mg}{l}$', color="blue")
    ax2.step(t, d, "blue", where='post', alpha=0.4)
    ax2.tick_params(axis='y', labelcolor="blue")
fig.tight_layout()  # otherwise the right y-label is slightly clipped
plt.xlabel('Time')
# plt.savefig('Protocol_realistic', dpi=300)
plt.show()
