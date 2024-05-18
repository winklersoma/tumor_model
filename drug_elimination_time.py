import numpy as np
import time


def exp_stoch(D0, r, time_step):
    """
  Simulates an exponential decay from D0 with rate r and returns with a list of lists consisting of
  - sampled cumulative random times
  - the states in the random times
  - the discrete times
  - the corresponding states (number of molecules)].
  """
    D = D0

    state = [D]
    time = [0]

    statistic = [D]
    time_stat = [0]
    while time_stat[-1] < 100 + time_step and D > 0:
        a = r * D
        Tau = np.random.exponential(scale=(1 / a))
        D -= 1
        state.append(D)
        time.append(time[-1] + Tau)

        disc_step = int((time[-1] - time_stat[-1]) / time_step)

        if disc_step > 0:
            for j in range(disc_step):
                time_stat.append(time_stat[-1] + time_step)
                statistic.append(D)

    return [time, state, time_stat, statistic]


D0 = 10000000
r = np.log(2) / 2
time_step = 0.05
start = time.time()
Data = []
repeat = 20
for i in range(repeat):  # The simulation will run 'repeat' number of times!
    d = exp_stoch(D0, r, time_step)
    Data.append(d)
end = time.time()
print('Average machine time spent on one loop:', (end - start) / repeat)
