import numpy as np
import matplotlib.pyplot as plt

# Constants
k = 8.617333262 * 10 ** -8  # Boltzmann constant [keV/K]
c = 29.9792458  # speed of light [cm/ns]
a = 8.563279012 * 10 ** 22  # radiation constant converted from [J/m^3-K^4] to [keV/cm^3-keV^4]


# Define variable problem parameters
sigma_a = 0.5  # absorption opacity [cm^-1]
T0_in_K = 300  # initial temperature in [K]
T0 = k * T0_in_K  # initial temperature in [keV]

# Zone length [cm]
x0 = 0
x1 = 1
Nx = 100  # Number of discrete positions to emit from
# Time step [ns]
t0 = 0.0
t1 = 1
Nt = 100
# Angles (even numbers only)
Nu = 100  # number of discrete angles

# From inputs, create position, time, and angle arrays for the deterministic method
particle_position = x0 + (np.arange(Nx) + 0.5) * (x1 - x0) / Nx
times_to_emit = t0 + (np.arange(Nt) + 0.5) * (t1-t0) / Nt
angles_to_emit = -1 + (np.arange(Nu) + 0.5) * 2 / Nu
print('Particle positions:', particle_position)
print('Times of emission:', times_to_emit)
print('Angles of emission:', angles_to_emit)

# Initial Calculations
total_particles_generated = Nx * Nt * Nu
total_energy_NxNu = a * T0 ** 4 * (x1-x0)  #
total_energy_NxNtNu = c * sigma_a * a * T0 ** 4 * (t1-t0) * (x1-x0)
e_particle = total_energy_NxNtNu / total_particles_generated  # initial energy per particle

# Store particle information for the deterministic method
particles = []

particle_counter = 1  # Counter for particle numbering
total_energy_deposited = 0  # Initialize total_energy_deposited
sum_remaining_energies_in_zone = 0  # Initialize sum_remaining_energies_in_zone

# Iterate through the different combinations for the deterministic method
for i in range(Nx):
    for t in times_to_emit:
        for angle in angles_to_emit:
            particle = {
                "Name": f"Particle{particle_counter}",
                "Initial_Position": particle_position[i],
                "Time_of_Emission": t,
                "Angle_of_Emission": angle
            }
            # Calculate the distance to census.
            s = c * (t1 - particle["Time_of_Emission"])

            # Calculate energy deposited
            energy_deposited = e_particle * np.exp(-1 * sigma_a * s)

            particle["Energy_Deposited"] = energy_deposited
            particle["Path Length"] = s

            particles.append(particle)
            particle_counter += 1

            # Add the energy deposited to total_energy_deposited
            total_energy_deposited += energy_deposited

            remaining_energy = e_particle - energy_deposited
            sum_remaining_energies_in_zone += remaining_energy


# Variables for the Monte Carlo method
particles_mc = []
particle_counter_mc = 1
total_energy_deposited_mc = 0
sum_remaining_energies_in_zone_mc = 0

# Iterate through particles for the Monte Carlo method
for _ in range(total_particles_generated):
    # Generate random position and time for the particle
    position = np.random.uniform(x0, x1)
    time = np.random.uniform(t0, t1)

    # Generate a non-zero random angle
    angle = 0
    while angle == 0:
        angle = np.random.uniform(-1, 1)

    # Calculate the distance to census.
    s = c * (t1 - time)

    # Calculate energy deposited
    energy_deposited = e_particle * np.exp(-1 * sigma_a * s)

    # Add the energy deposited to total_energy_deposited
    total_energy_deposited_mc += energy_deposited

    remaining_energy = e_particle - energy_deposited
    sum_remaining_energies_in_zone_mc += remaining_energy

    # Store particle information for the Monte Carlo method
    particle_mc = {
        "Name": f"Particle{particle_counter_mc}_mc",
        "Initial_Position": position,
        "Time_of_Emission": time,
        "Angle_of_Emission": angle,
        "Energy_Deposited": energy_deposited,
        "Path Length": s
    }

    particles_mc.append(particle_mc)
    particle_counter_mc += 1

# Print Problem Parameters

print('****************************')
print('**** Problem Parameters ****')
print('****************************')
print()
print(' Constants used:'
      '\n Boltzmann constant, k = {} [keV/K]'
      '\n speed of light, c = {} [cm/ns]'
      '\n radiation constant, a = {} [keV/cm^3-keV^4]'
      '\n'
      '\nOther Properties:'
      '\nabsorption opacity, sigma_a = {} [cm^-1]'
      ''.format(k, c, a, sigma_a)
      )

print('The zone where particles are sourced is from {} to {} cm.'.format(x0, x1))
print(f'The initial temperature in the zone is {T0_in_K} K.')
print('{} particles were generated from {} positions in the source zone,'
      ' {} time discretizations,'
      ' and {} discrete angles'
      ' in a {} ns time step.'.format(total_particles_generated, Nx, Nt, Nu, t1-t0))

# Print results for deterministic method
print()
print('Results for Deterministic Sourcing Method')
print()
print('Total source particles generated:', total_particles_generated)
print('Total Energy of initial particles [NxNu] = a * T0^4 * Vz =', total_energy_NxNu)
print('Total Energy of all source particles [NxNuNt] = c * sigma_a * a * T0^4 * delta_t * Vz =', total_energy_NxNtNu)
print('Energy deposited in the material =', total_energy_deposited)
print('Energy of radiation remaining =', total_energy_NxNtNu - total_energy_deposited)
print('Ideal Energy deposited = c * sigma_a * a * T0^4 * delta_t * Vz=', total_energy_NxNtNu)
print('Ideal Energy of particles remaining = a * T0^4 * Vz=', total_energy_NxNu)


# Print results for Monte Carlo Sourcing Method
print()
print('Results for Monte Carlo Sourcing Method')
print()
print('Total source particles generated:', total_particles_generated)
print('Total Energy Deposited [keV] for Monte Carlo Method =', total_energy_deposited_mc)
print('Remaining radiation energy [keV] for Monte Carlo Method =', sum_remaining_energies_in_zone_mc)
print('Final Energy [keV] for Monte Carlo Method =', total_energy_deposited_mc + sum_remaining_energies_in_zone_mc)
energy_difference_mc = total_energy_NxNtNu - (total_energy_deposited_mc + sum_remaining_energies_in_zone_mc)
print('Energy Difference for Monte Carlo Method =', energy_difference_mc)
print('a*T^4* Vz= ', a * T0 ** 4 * (x1-x0))
