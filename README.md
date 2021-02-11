# Ising-model-for-Triangular-lattice
MATLAB program to esimate average Curie or critical Temperature (Tc) and uncertainty in the value of estimated Tc (in Kelvin) for a given DMS with a triangular lattice. 

Background: In the code, DMS is a hexagonal molybdenum disulfide monolayer (MoS2), a semiconductor material, with a small concentration of magnetic impurities of manganese (Mn) dopants. In MoS2, Mo atoms assume a triangular lattice, and to get DMS some Mo atoms are randomly replaced or doped with Mn dopants. 

DMS is modeled as a triangular lattice with periodic boundary conditions (PBCs) in two in-plane directions. Each site in the lattice represents Mo or Mn atom with spin (s). Within the 2D Ising model, Mo atoms have fixed spin, s = 0, and magnetic Mn atoms have variable spins, s = +1 or -1.

Within the 2D Ising model, in absence of an external magnetic field, the total energy of the system includes interaction terms between neighboring spins. To keep calculation tractable, long-range spin-spin interactions between atoms within the first four nearest neighbor shells are considered. Code implements the Metropolis Monte Carlo scheme (MMC), to minimize the total energy of the system at temperature, T, by explicitly flipping the spin of randomly chosen Mn atom. 

At a certain T, code calculates four quantities: total energy, specific heat, specific magnetization, and magnetic susceptibility. To include the effect of configurational disorder of Mn atoms in a lattice, code repeats a calculation multiple times with different random initial distribution of Mn atoms and then calculates the average of all thermodynamic/magnetic properties over these distinct runs. Using this scheme, code increases temperature, T, from 0.05 to 0.70 in the steps of 0.01 (T is in non-dimensional units) and stores average thermodynamic/magnetic properties in a list at each T. 

In the absence of an external magnetic field, the system undergoes a phase transition, from ferromagnetic phase to paramagnetic phase, when the magnetic susceptibility diverges at the Tc. Code finds the Tc and outputs average Tc and uncertainty in the value of estimated Tc in Kelvin.

Code usage:
The main function final_ising () calculates the Tc of DMS with a triangular lattice, and to do so it needs three arguments: 
(1)	Size of the lattice (L). Example: L=20 will generate 20 x 20 supercell of a triangular lattice.
(2)	The fraction of magnetic dopant in the lattice (conc_Mn). Example: conc_Mn = 0.10 corresponds to 10% Mn dopants in lattice.
(3)	Distinct runs (repeat). Example: if repeat=10, code will perform 10 calculations with a distinct random initial distribution of Mn atoms in lattice and calculate average thermodynamic/magnetic properties over these runs. 
            
The final_ising () initializes parameters and variables of the 2D Ising model in non-dimensional units and then implements the Metropolis Monte Carlo scheme to minimize the total energy of the system and estimate the average Tc and uncertainty in Tc. Parameters and interactions terms needed for energy calculation comes from “Physical Review B 87.19 (2013): 195201”.

The output of final_ising (): Average Tc and uncertainty in Tc in Kelvin (K). 

Within the final_ising (), two functions f1_shell () and f2_shell () makes a list of 1st and 2nd nearest neighbors, respectively, of a chosen atom/site in a triangular lattice with PBCs in two in-plane directions. f2_shell () is further used to find 3rd and 4th nearest neighbors in successive shells. 

f1_shell () takes three arguments: 
(1)	x and y coordinates of the chosen atom in a triangular lattice (site)
(2)	List, S, which contains the spin of each atom
(3)	Size of the lattice (L)

The output of f1_shell (): List of coordinates of six 1st nearest neighbors of the chosen atom in a triangular lattice. 

The function f2_shell () uses f1_shell () to find six 1st nearest neighbors of each atom provided in first_next_shell_nn, and it then removes repeated entries, to finally output the list of unique 2nd nearest neighbors.

f2_shell () takes four arguments: 
(1)	x and y coordinates of the chosen atom in a triangular lattice (site)
(2)	A list (first_next_shell_nn) contains x and y coordinates of atoms  
(3)	List, S, which contains the spin of each atom
(4)	Size of the lattice (L)

The output of f2_shell (): List of 2nd nearest neighbors of the chosen atom in a triangular lattice
