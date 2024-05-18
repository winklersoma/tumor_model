from scipy.optimize import minimize

import functions as fun
from smoothing import *

np.savetxt('test.txt', np.array([0, 1]), fmt='%.1e', delimiter=',')

vita_R = s_c8161_R[:40]
vita_Y = s_c8161_Y[:40]
vita_G = s_c8161_G[:40]


def fn_min(params):
    t, r, y, g, m = fun.cc_3_logdelay(params)
    simu_R, simu_Y, simu_G = m + r, y, g
    return np.sum((simu_R - vita_R) ** 2) + np.sum((simu_Y - vita_Y) ** 2) + np.sum((simu_G - vita_G) ** 2)  # *


parameters = np.array([4.187304626818229e+01,
                       6.409128978848699e+00,
                       1.602263973162360e+00,
                       1.174141163898582e+01,
                       8.133619971889119e-01,
                       1.236030451168313e+01,
                       1.138172769181080e+00,
                       9.714652235883265e-01,
                       9.886208021886280e-01,
                       6.688432723645765e-01,
                       8.242534563840548e-01,
                       1.337292772552692e+00,
                       4.399728563416122e-01])

res = minimize(fn_min, parameters, method='Nelder-Mead',
               options={'maxiter': 500})

# bounds=[[0.1,51],[3,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2],[5,15],[0.1,2]]

print(res)
print(res.x)

np.savetxt('results2.txt', np.array(res.x), fmt='%.15e', delimiter=',')
