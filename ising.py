import sys
import numpy as np
import matplotlib.pyplot as plot
from scipy.constants import k

# params for simulation
runs = 9**7
frames = 500

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

def mcmove(beta):
    rejections = 0

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

            if i % frames == 0:
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
    print(f'{rejections} rejections')

magvector = []
N = rows*cols

for T in Ts:
    lattice = np.random.choice([-1, 1], size=(rows, cols))
    mcmove(1/T)

    netmag = np.sum(lattice) / N 

    print(f'Magnetization = {netmag}')
    magvector.append(abs(netmag))

plot.scatter(Ts, magvector)
plot.xlabel('Temperature J/K')
plot.ylabel('Magnetization')
plot.title('Temp v Mag for 2D ising model')
plot.show()
