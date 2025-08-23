import numpy as np
import matplotlib.pyplot as plot
from scipy.constants import k

# params for simulation
runs = 10**7
frames = 500

# params for physical system
j = 1 #align, -1 for anti-align
T = 1 #room temp J/k_b
# k_b = k # Boltzmann constant: for some reason this results in z values below being only 0.0 or 1.0 (possible float arithmetic issue?)
k_b = 1 #fix

# lattice params
rows = 200
cols = 200
lattice = np.random.choice([-1, 1], size=(rows, cols))

#consts
beta = 1/(k_b*T)

def hamiltonian(latt, sigy, sigx):
    sumoverlatt = 0

    if 0 <= sigy-1: sumoverlatt += latt[sigy-1][sigx] * latt[sigy][sigx]
    if sigy+1 < cols:  sumoverlatt += latt[sigy+1][sigx] * latt[sigy][sigx] 
    if 0 <= sigx-1: sumoverlatt += latt[sigy][sigx-1] * latt[sigy][sigx]
    if sigx+1 < rows:  sumoverlatt += latt[sigy][sigx+1] * latt[sigy][sigx]
    return -j*sumoverlatt

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
            axes.text(0.5, -0.05, f"Step: {i}",
                transform=axes.transAxes,
                ha='center', va='top', fontsize=12)
            plot.draw()
            plot.pause(0.01)
except KeyboardInterrupt:
    print("Simulation stopping...")
finally:
    plot.ioff()
    plot.show()

plot.ioff()
plot.show()
print(f'{rejections} rejections')

