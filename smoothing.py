import matplotlib.pyplot as plt

from data import *

weights = np.array([1 / 9, 2 / 9, 3 / 9, 2 / 9, 1 / 9])

s_wm983c_R = [wm983c[0][:2]]
s_wm983c_Y = [wm983c[1][:2]]
s_wm983c_G = [wm983c[2][:2]]
length_wm = len(wm983c[0]) - 4

s_lu1205_R = [lu1205[0][:2]]
s_lu1205_Y = [lu1205[1][:2]]
s_lu1205_G = [lu1205[2][:2]]
length_lu = len(lu1205[0]) - 4

s_c8161_R = [c8161[0][:2]]
s_c8161_Y = [c8161[1][:2]]
s_c8161_G = [c8161[2][:2]]
length_c8 = len(c8161[0]) - 4

for i in range(0, length_wm):
    s_wm983c_R = np.append(s_wm983c_R, np.sum(weights * wm983c[0][i:i + 5]))
    s_wm983c_Y = np.append(s_wm983c_Y, np.sum(weights * wm983c[1][i:i + 5]))
    s_wm983c_G = np.append(s_wm983c_G, np.sum(weights * wm983c[2][i:i + 5]))

s_wm983c_R = np.append(s_wm983c_R, [wm983c[0][-2:]])
s_wm983c_Y = np.append(s_wm983c_Y, [wm983c[1][-2:]])
s_wm983c_G = np.append(s_wm983c_G, [wm983c[2][-2:]])

for i in range(0, length_lu):
    s_lu1205_R = np.append(s_lu1205_R, np.sum(weights * lu1205[0][i:i + 5]))
    s_lu1205_Y = np.append(s_lu1205_Y, np.sum(weights * lu1205[1][i:i + 5]))
    s_lu1205_G = np.append(s_lu1205_G, np.sum(weights * lu1205[2][i:i + 5]))

s_lu1205_R = np.append(s_lu1205_R, [lu1205[0][-2:]])
s_lu1205_Y = np.append(s_lu1205_Y, [lu1205[1][-2:]])
s_lu1205_G = np.append(s_lu1205_G, [lu1205[2][-2:]])

for i in range(0, length_c8):
    s_c8161_R = np.append(s_c8161_R, np.sum(weights * c8161[0][i:i + 5]))
    s_c8161_Y = np.append(s_c8161_Y, np.sum(weights * c8161[1][i:i + 5]))
    s_c8161_G = np.append(s_c8161_G, np.sum(weights * c8161[2][i:i + 5]))

s_c8161_R = np.append(s_c8161_R, [c8161[0][-2:]])
s_c8161_Y = np.append(s_c8161_Y, [c8161[1][-2:]])
s_c8161_G = np.append(s_c8161_G, [c8161[2][-2:]])

###########################    Visualization    #################################

plt.scatter(np.arange(0, len(wm983c[0]) * 0.25, 0.25), wm983c[0], alpha=0.2, c="red")
plt.scatter(np.arange(0, len(wm983c[1]) * 0.25, 0.25), wm983c[1], alpha=0.2, c="orange")
plt.scatter(np.arange(0, len(wm983c[2]) * 0.25, 0.25), wm983c[2], alpha=0.2, c="green")

plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), s_wm983c_R, linewidth=2.0, alpha=1, c="red")
plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), s_wm983c_Y, linewidth=2.0, alpha=1, c="orange")
plt.plot(np.arange(0, len(wm983c[2]) * 0.25, 0.25), s_wm983c_G, linewidth=2.0, alpha=1, c="green")
plt.title('The R, Y, G curves after the smoothing process')
plt.xlabel('Time')
plt.ylabel('Number of cells')
plt.savefig('smoothing.png', dpi=300)
plt.show()
