# 2D Magnetic Phase Transition Simulation Ising Model, using Monte Carlo
### Theory
#### 1. 1D Ising Model
    * We first calculate the Hamiltonian: H = -J * sum-over-lattice(sigma(i) * sigma(j)), where
    sum-over-lattice evaluates the equation N times. Therefore, there are a possible 2^N 
    possible spin configurations. 
    * We then need to weigh Boltzmann factors. The probability of a configuration appearing is not
    uniform, its weighted by e^-BH, where H is the Hamiltonian we calculated and B=1/(Kb * T)
    * Now, we evaluate the partition function: Z = sum-over-spinstates(e^-BH(spinstate))
    * NOTE: In a 1D Ising Model, there is no phase transition. However, applying these same equations to
    our 2D Ising Model shows that there is a phase transition (change in magnetization without 
    any external field).
#### 2. 2D Ising Model
    * We could prove that a 2D Ising Model undergoes a phase transition above a certain critical 
    temperature using analytical methods (see Lans Onsanger's Transfer Matrix method), but we 
    are more interested in using a Monte-Carlo algorithm to numerically simulate this 2D ising
    model system.

### Model
    1. We initialize a lattice with randomized spins of -1 and 1. P(s -> s') describes the probability 
    of a spin state to swap; i.e., go from -1 to 1 or vice versa.
    2. We describe the probability function as:
        P(s -> s')  = 1                         ... if H(s') < H(s)
                    = e^(-B * (H(s')-H(s)))     ... otherwise
    3. We create a nice graphical representation of spin states that evolves with time to (hopefully) 
    observe a phase transition
    * NOTE: We don't actually need to compute the hamiltonian for the entire lattice; once we roll a
    random cell, we can compute Hi before the flip and Hf after, and compare those to get the exact
    same result for delta-H
