import sys
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plt
import signal
from numba import njit

# Params
sweeps_per_cell = 5 * (10**4)
frames = 1500
samples_per_temp = 10

gui = True

# Physical params
j = 1
T_min = 1.0
T_max = 3.5
T_n = 50

rows = 50
cols = 50

Ts = np.linspace(T_min, T_max, T_n)

# Physics
def hamiltonian(latt, y, x):
    return -j * (
        latt[(y+1) % rows][x] * latt[y][x] +
        latt[(y-1) % rows][x] * latt[y][x] +
        latt[y][(x+1) % cols] * latt[y][x] +
        latt[y][(x-1) % cols] * latt[y][x]
    )


def mcmove_gui(beta, lattice):
    runs = sweeps_per_cell * rows * cols;

    plt.ion()
    fig, ax = plt.subplots()
    img = ax.imshow(lattice, cmap="spring", interpolation="nearest")
    ax.axis("off")

    for i in range(runs):

        x = np.random.randint(0, cols)
        y = np.random.randint(0, rows)

        hi = hamiltonian(lattice, y, x)
        dH = -2 * hi

        if dH <= 0 or np.random.rand() < np.exp(-beta * dH):
            lattice[y][x] *= -1

        if (i % frames == 0):
            img.set_data(lattice)
            ax.set_title(f"T = {1/beta:.2f}, step = {i:,}")
            fig.canvas.draw_idle()
            plt.pause(0.001)

    plt.ioff()
    img.set_data(lattice)
    ax.set_title(f"Final state, T = {1/beta:.2f}")
    plt.show()

@njit(cache=True)
def mcmove_nogui(beta, lattice, j):
    runs = sweeps_per_cell * rows * cols;
    for i in range(runs):
        x = np.random.randint(0, cols)
        y = np.random.randint(0, rows)
        hi = -j * (
            lattice[(y+1) % rows, x] * lattice[y, x] +
            lattice[(y-1) % rows, x] * lattice[y, x] +
            lattice[y, (x+1) % cols] * lattice[y, x] +
            lattice[y, (x-1) % cols] * lattice[y, x]
        )
        dH = -2 * hi
        if dH <= 0 or np.random.rand() < np.exp(-beta * dH):
            lattice[y, x] *= -1


def simulatefortemp(T, show):
    try:
        beta = 1 / T

        mag = []
        for i in range(samples_per_temp):
            lattice = np.random.choice([-1, 1], size=(rows, cols)).astype(np.int32)
            if show:
                mcmove_gui(beta, lattice)
            else:
                mcmove_nogui(beta, lattice, j)

            mag.append(abs(np.sum(lattice)) / (rows * cols))
        print(f'Temperature: {T} | Magnetization: {np.mean(mag)}\n')
        return T, mag

    except KeyboardInterrupt:
        return None

def worker_init():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

if __name__ == "__main__":

    mp.freeze_support()

    nproc = mp.cpu_count()
    chunks = [Ts[i:i+nproc] for i in range(0, len(Ts), nproc)]

    results = []

    for chunk in chunks:
        args = [(chunk[0], gui)] + [(T, False) for T in chunk[1:]]
        with mp.Pool(processes=nproc, initializer=worker_init) as pool:
            results_async = pool.starmap_async(simulatefortemp, args)
            try:
                while not results_async.ready():
                    results_async.wait(timeout=0.5)
                res = results_async.get()
                results.extend([r for r in res if r is not None])
            except KeyboardInterrupt:
                pool.terminate()
                pool.join()
                print("\nStopped by user")
                sys.exit()

    # -----------------------
    # Plot results
    # -----------------------
    results.sort(key=lambda x: x[0])

    Ts_sorted = [r[0] for r in results]
    mags = [r[1] for r in results]

    mags_mean = [np.mean(r[1]) for r in results]
    mags_std  = [np.std(r[1])  for r in results]

    plt.figure()
    plt.errorbar(Ts_sorted, mags_mean, yerr=mags_std,
                 fmt='o', capsize=4, capthick=1, elinewidth=1, markersize=5)
    plt.xlabel("Temperature")
    plt.ylabel("Magnetization")
    plt.title("2D Ising Model: Magnetization vs Temperature")
    plt.show()