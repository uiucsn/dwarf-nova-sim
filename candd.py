import numpy as np
from matplotlib import pyplot as plt

def data_loading(name):
    dff = np.genfromtxt(name, names='t,m,err')
    t = dff['t']
    m = dff['m']
    return t, m

def find_cand(t, m, mean, sigma, n, n_below, t_before, t_aft):
    diff_bool_idx = np.diff(np.asarray(m >= mean - n * sigma, dtype=int))
    left_idx = np.where(diff_bool_idx == -1)[0]
    right_idx = np.where(diff_bool_idx == 1)[0]

    candidates = []
    for i_seq, (left, right) in enumerate(zip(left_idx, right_idx)):
        print(i_seq)
        print((left, right))
        sl = slice(left + 1, right + 1)
        print(sl)

        magn = m[sl]

        #     magnerr = err[sl]
        time = t[sl]

        duration = time[-1] - time[0]

        if duration >= 20 or duration < 4:
            continue

        for i in range(1, left + 1):
        # for i in range(left - 1, -1, -1):
            # t[i]

            if t[left] - t[left - i] > t_before:
                i += 1
                break
            if m[left - i] > mean + n_below * sigma:
                break
            if m[left - i] < m[left - i + 1]:
                i += 1
                break
        time = np.insert(time, 0, t[left - i:left])
        magn = np.insert(magn, 0, m[left - i:left])
        # points after outburst
        # TODO: rewrite as left
        for j in range(1, len(time) - right):
        # for j in range(right + 1, len(time)):
            # j = 1 .. right - 1
            # right + j = right + 1 .. 2 * right - 1

            if t[right + j] - t[right] > t_aft:
                break
            elif (m[right + j] <= mean + n_below * sigma) and (m[right + j] >= m[right + j - 1]):

                time = np.append(time, t[right + j])
                magn = np.append(magn, m[right + j])

        candidates.append((i_seq, time, magn))
    return candidates

def plot_cand(t, m, candidates, mean, sigma, n, n_below):
    fig, axes = plt.subplots()
    axes.axhline(y=mean)
    axes.axhline(y=mean - n * sigma, c='g', linestyle='dotted')
    axes.axhline(y=mean + n_below * sigma, c='g', linestyle='dotted')
    # print all:
    axes.scatter(t, m, c='c', marker="x")

    for i in range(len(candidates)):
        t1 = candidates[i][1]
        m1 = candidates[i][2]
        axes.scatter(t1, m1, c='r', marker="x")

    axes.invert_yaxis()
    plt.show()

def main():
    t, m = data_loading("phot/OGLE-BLG-DN-0002.dat")
    sigma = np.std(m)
    mean = np.mean(m)
    n = 2
    n_below = 0.5
    t_before = 3
    t_aft = 6
    candidates = find_cand(t, m, mean, sigma, n, n_below, t_before, t_aft)
    plot_cand(t, m, candidates, mean, sigma, n, n_below)


if __name__ == '__main__':
    main()
