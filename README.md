# 2D Magnetic Phase Transition Simulation Ising Model, using Monte Carlo
### Theory
#### 1. 1D Ising Model
- We first calculate the Hamiltonian: $H=-J \cdot \sum_{<i,j>}\sigma_i\sigma_j$, where
the $\sum_{<i,j>}$ term evaluates the equation $N$ times. Therefore, there are a possible $2^N$ 
possible spin configurations. 
- We then need to weigh Boltzmann factors. The probability of a configuration appearing is not
uniform, its weighted by $e^{- \beta H}$, where $H$ is the Hamiltonian we calculated and $$\beta = \frac{1}{k_B \cdot T}$$ 
- Now, we evaluate the partition function: $Z = \sum_{s_i}e^{- \beta H(s_i)}$, where $s_i$ encodes all
possible configurations of the lattice.
- NOTE: In a 1D Ising Model, there is no phase transition. However, applying these same equations to
our 2D Ising Model shows that there is a phase transition (change in magnetization without 
any external field).
#### 2. 2D Ising Model
- We could prove that a 2D Ising Model undergoes a phase transition above a certain critical 
temperature using analytical methods (see [Transfer Matrix method](https://en.wikipedia.org/wiki/Ising_model#Two_dimensions)), but we 
are more interested in using a Monte-Carlo algorithm to numerically simulate this 2D ising
model system.

### Model
#### 1. Procedure
1. We initialize a lattice with randomized spins of -1 and 1. $P(\sigma\to\sigma')$ describes the probability 
of a spin state to swap; i.e., go from -1 to 1 or vice versa.
2. We describe the probability function as:

$$P(\sigma \to \sigma') = \begin{cases} 
1, & \text{if} \\ H(\sigma) < H(\sigma') \\ 
e^{- \beta \Delta H}, & \text{otherwise} \end{cases}$$

3. We create a nice graphical representation of spin states that evolves with time to (hopefully) 
observe a phase transition
- NOTE: We don't actually need to compute the hamiltonian for the entire lattice; once we roll a
random cell, we can compute $H_i$ before the flip and $H_f$ after, and compare those to get the exact
same result for $\Delta H$.
#### 2. Parameters
- Simulation parameters:
  * `runs`: Number of Monte Carlo runs, (default) `5*10**7`
  * `gui`: Enable GUI to view lattice evolution. Because of multiproc constraints, we only show Main thread output. (default) `False`
  * `frames`: Number of frames to update the simulation, only if GUI is enabled. `500`
- Physical parameters:
  * `J`: Alignment operator; 1 for alignment, -1 for anti-alignment. (default) `1`
  * `T_min`: Start temperature for simulation. (default) `0.5`
  * `T_max`: End temperature for simulation. (default) `4.5`
  * `T_n`: # of simulation steps. (default) `20`
