import numpy as np
from matplotlib import pyplot as plt
import os

def data_loading(name):
    dff = np.genfromtxt(name, names='t,m,err')
    t = dff['t']
    m = dff['m']
    return t, m

def find_cand(t, m, mean, sigma, n, n_below, t_before, t_aft):
    diff_bool_idx = np.diff(np.asarray(m >= mean - n * sigma, dtype=int))
    left_idx = np.where(diff_bool_idx == -1)[0]
    right_idx = np.where(diff_bool_idx == 1)[0]
    i = None
    j = None
    candidates = []
    for i_seq, (left, right) in enumerate(zip(left_idx, right_idx)):

        sl = slice(left + 1, right + 1)


        magn = m[sl]


        time = t[sl]

        if time.size < 1:
            continue

        duration = time[-1] - time[0]

        if duration >= 20 or duration < 4:
            continue

        for i in range(left - 1, -1, -1):

            if t[left] - t[i] > t_before:
                i += 1
                break


            if m[i] > mean + n_below * sigma:
                break
            # if m[left] < m[i+1]:

            #     i += 1
            #     break

        time = np.insert(time, 0, t[i:left])
        magn = np.insert(magn, 0, m[i:left])
        # points after outburst
        #

        for j in range(len(time) - right):
            # j = 1 .. right - 1
            # right + j = right + 1 .. 2 * right - 1


            if t[right + j] - t[right] > t_aft:
                break
            elif (m[right + j] <= mean + n_below * sigma) :

                time = np.append(time, t[right + j])
                magn = np.append(magn, m[right + j])
        
        if time.size < 8:
            continue
        
        if time[-1]-time[0]<15:
            continue
        candidates.append((i_seq, time, magn))
    return candidates

def plot_cand(t, m, candidates, mean, sigma, n, n_below, name):
    fig, axes = plt.subplots()
    axes.axhline(y=mean)
    axes.axhline(y=mean - n * sigma, c='g', linestyle='dotted')
    axes.axhline(y=mean + n_below * sigma, c='g', linestyle='dotted')

    axes.scatter(t, m, c='c', marker="x")
    # axes.set_xlim(6300, 6400)
    for i in range(len(candidates)):
        t1 = candidates[i][1]
        m1 = candidates[i][2]
        axes.scatter(t1, m1, c='r', marker="x")

    axes.invert_yaxis()
    plt.xlabel('time /mjd')
    plt.ylabel('magnitude')

    basena = os.path.basename(name)
    title, extension = os.path.splitext(basena)
    plt.title(title)
    direc = 'pictures'
    os.makedirs(direc, exist_ok=True)
    plt.savefig(os.path.join(direc, f'{title}.png'))
    plt.close()
    # plt.show()
def main():
    for i in range(1, 50):
        # i = 11
        numb =  "%04d" % i
        name = "phot/OGLE-BLG-DN-" + str(numb) + ".dat"

        t, m = data_loading(name)
        sigma = np.std(m)
        mean = np.mean(m)
        n = 2
        n_below = 0.5
        t_before = 3
        t_aft = 6
        candidates = find_cand(t, m, mean, sigma, n, n_below, t_before, t_aft)
        plot_cand(t, m, candidates, mean, sigma, n, n_below, name)



if __name__ == '__main__':
    main()
