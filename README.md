# Fluid System

The main code is in the FluidGUI folder.

## Summary

For my project, I decided to simulate a 3D fluid system using the SPH framework. Within this system, users
can control physical forces, tank dimensions, and various particle attributes on a container of fluids.
The system has been optimized using a Spatial Hashing algorithm for fast neighbour searching of particles
in a neighbourhood of one another.

## Main Framework

+ ![image missing](images/Main_Class_Setup.PNG)

### This system is comprised of multiple parts which are executed before and after the simulation:
#### On initialization:
- Hashing particles, to find particle neighbours of a given particle in the simulation
- adding these particles based on their hash key to individual cells in a partitioned grid of particles
- calculating the forces on each particle based on their neighbour's attributes
- setting the sum of forces, acceleration, velocity and position of each particle
- correcting velocity using XSPH velocity correction on each particle
- checking for boundary/tank collisions per particle and adjusting positions and velocities based on these

#### While simulation runs
- The hash table is dynamically updated based on updated particle positions
- Reassigning existing particles to their appropriate cells using an update hashing method
- Retrieving neighbour id's from particles stored in the same cell
- Recomputing the sum of forces, acceleration, velocities and positions of each particle
- restart this process per frame.

## Spatial Hashing Algorithm

- For faster performance, I used spatial hashing to set unique hash keys to particles based on their positions,
  and whenever a particle was found, I set it's hash key and added its id to an unordered map. This map stored
  all the id's of particles with the same hash key into a unique cell. Whenever I needed to access neighbour 
  positions from a given cell, I created a function to get particle id's from a given cell in the hash table.
- For larger scale simulations, this can give more accurate results and performs at O(n) complexity as opposed
  to the brute force approach to finding particle neighbours.

+ ![image missing](images/hashing_algorithm.PNG)

+ ![image missing](images/hashing_ids.PNG)

+ ![image missing](images/next_prime_calc.PNG)


## SPH Framework
### Class Framework

The entire project is split into a Tank, System, Particle, NGLScene and MainWindow class:


+ ![image missing](images/Main_Class_Variables.PNG)

#### Tank Class:
In the tank class, I specify getters and setters to act upon the fluid tank holding the particles in the 
scene. These include operations to set tank dimensions and positions in the scene space.

#### System Class
The system class performs all the operations on an AbstractVAO object of GL_POINTS, which stores all particle 
data inherited from the Particle class into an std::vector and computes forces and updates particle positions.
This class also inherits slots from the NGLScene to update particle and system attributes from user inputs. 


##### Force Calculations
+ ![image missing](images/physical_parameter_defaults.PNG)

+ ![image missing](images/pressure_force_calculation.PNG)

+ ![image missing](images/mass_density_calculation.PNG)

+ ![image missing](images/surface_curvature_calculation.PNG)

+ ![image missing](images/viscosity_force_calculation.PNG)

+ ![image missing](images/position_velocity_calculation.PNG)

#### Particle Class
The particle class is a placeholder for each particle including getters and setters for particle id's, hashkeys
and individual std::vector neighbour lists storing the id's of particle neighbours to be computed later on by
the System class.

#### NGLScene Class
This class is used to visualize the system in an NGLScene and actually contains the operations to draw and update
the particle positions per frame. Separate ...func methods exist in this class to send user data from the main 
window class to the System class.

#### Main Window Class
The main window class is where the Signals and Slots are first created to bind user data to the NGLScene, which
subsequently modifies System and Particle data in the inheriting classes.


### Class Implementation

+ ![image missing](images/System_Class_Setup.PNG)

To implement this system, I was inspired by a master's paper that tackled a fluid simulation in NGL, but with some
distinct differences in the class setup. Instead of calculating per particle forces in an independent Particle class, 
I used the Particle class as a placeholder for per particle attributes. All the computations were instead performed
through the System class which inherited from the Particle class.

On update and initialization calls in the system class, the aforementioned simulation steps are performed per particle
once per frame in a for loop during update and initialization of the class. Each force is calculated with a single id
reference parameter, which retrieves particle data from its id, and calculates individual and neighbour forces acting
upon that particle from its id.

The pressure force and mass density functions also include an update pressure and density call before calculations are
performed to compute the pressure force on a particle.

Finally, a sum forces function exists, which updates all the forces on that particle from its id, which is then used to
calculate its position, velocity and acceleration.

## Kernel Weighting Functions

A crucial aspect of this implementation is the Kernel functions. Each kernel function type also has a gradient and 
laplacian operator depending on the type of force to be calculated. Different kernel types exist in this simulation, not
least the Poly 6 Kernel, Spiky Kernel and a viscosity kernel, adapted for accurate viscosity calculations.

The kernel approximations are used as weighting factors when determining forces per particle based on their 
neighbour's attributes.

For code reusability, I created 3 functions for kernel, kernel gradient and kernel laplacian calculations per particle.
Within the kernel function, a switch statement is used to vary between the poly 6, spiky and viscosity kernel based 
on the neccessary kernel type for the calculation. The same logic applied in both the kernel gradient and laplacian 
functions.

It was important to normalize values from the kernel function, as the Weighing kernel is said to evaluate between 0
and 1, and in my first test trials of the simulation, the values from the kernels would explode creating numerical instabilities
in the calculations.

## Final Video Link
[![Fluid Sim Link here.](https://img.youtube.com/vi/watch?v=6YAEiR-Dszc/0.jpg)](https://www.youtube.com/watch?v=6YAEiR-Dszc)

### References
BlankNGL - Used as a starting point for my system

Fluid Master's Paper - https://nccastaff.bournemouth.ac.uk/jmacey/MastersProject/MSc16/15/thesis.pdf

Mueller Paper - https://matthias-research.github.io/pages/publications/sca03.pdf


