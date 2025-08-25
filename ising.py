import sys
import numpy as np
import multiprocessing as mp
import matplotlib.pyplot as plot
from scipy.constants import k

# params for simulation
runs =  5 * (10**7)
frames = 500
gui = False

# params for physical system
j = 1 #align, -1 for anti-align
T_min = 0.5 
T_max = 4.5
T_n = 20
# k_b = k # Boltzmann constant: for some reason this results in z values below being only 0.0 or 1.0 (possible float arithmetic issue?)
k_b = 1 #fix

# lattice params
rows = 100
cols = 100

#consts
Ts = np.linspace(T_min, T_max, T_n)

def hamiltonian(latt, sigy, sigx):
    sumoverlatt = 0

    sumoverlatt += latt[(sigy+1) % rows][sigx] * latt[sigy][sigx] + latt[sigy][(sigx+1) % cols] * latt[sigy][sigx]
    sumoverlatt += latt[(sigy-1) % rows][sigx] * latt[sigy][sigx] + latt[sigy][(sigx-1) % cols] * latt[sigy][sigx]
    return -j*sumoverlatt

def mcmove(beta, plotprog):
    plot.ion()
    figure, axes = plot.subplots()

    try: 
        for i in range(runs):
            randx, randy = np.random.randint(0,rows), np.random.randint(0,cols)

            hi = hamiltonian(lattice, randy, randx)
            dH = -2*hi

            if dH <= 0:
                lattice[randy][randx] *= -1
            else:
                pi_s = np.random.rand()
                z = np.exp(-beta * dH)

                if pi_s < z:
                    lattice[randy][randx] *= -1

            if plotprog and i % frames == 0:
                axes.clear()
                axes.imshow(lattice, cmap="spring", interpolation='nearest')
                axes.axis("off")
                axes.text(0.5, -0.05, f"Step: {i}, T: {1/beta}",
                    transform=axes.transAxes,
                    ha='center', va='top', fontsize=12)
                plot.draw()
                plot.pause(0.01)
    except KeyboardInterrupt:
        print("Simulation stopping...")
        sys.exit()

    plot.ioff()
    plot.close()

def simulatefortemp(T, plotprog):
    N = rows*cols
    global lattice
    lattice = np.random.choice([-1, 1], size=(rows, cols))
    mcmove(1/T, plotprog)
    netmag = abs(np.sum(lattice) / N) 
    print(f'Magnetization = {netmag} for Temp = {T}')
    return (T, netmag)


nproc = mp.cpu_count()
chunks = [Ts[i:i+nproc] for i in range(0, len(Ts), nproc)]
results = []

for chunk in chunks:
    args = [(chunk[0], gui)] + [(T, False) for T in chunk[1:]]
    with mp.Pool(processes=mp.cpu_count()) as pool:
        chresults = pool.starmap(simulatefortemp, args)
    results.extend(chresults)

if __name__ == '__main__':
    results.sort(key= lambda x: x[0])
    Ts_sort = [r[0] for r in results]
    magvector = [r[1] for r in results]

    plot.scatter(Ts_sort, magvector)
    plot.xlabel('Temperature J/K')
    plot.ylabel('Magnetization')
    plot.title('Temp v Mag for 2D ising model')
    plot.show()
