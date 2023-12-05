import matplotlib.pyplot as plt
import numpy as np

# Constants
k = 8.617333262e-8  # Boltzmann constant [keV/K]
c = 29.9792458      # speed of light [cm/ns]
a = 8.563279012e22  # radiation constant converted from [J/m^3-K^4] to [keV/cm^3-keV^4]

# Define variable problem parameters
sigma_a = 0.5     # absorption opacity [cm^-1]
T0_in_K = 300     # initial temperature [K]
T0 = k * T0_in_K  # initial temperature [keV]

# Source zone length [cm]
x0 = 0
x1 = 1
Vz = x1 - x0  # zone volume
# Time step [ns]
t0 = 0
t1 = 0.03


def simulate_deterministic(num_particles):
    Nx = num_particles
    Nt = num_particles
    Nu = num_particles
    # Create position, time, and angle arrays for the deterministic method
    particle_position = x0 + (np.arange(Nx) + 0.5) * Vz / Nx
    times_to_emit = t0 + (np.arange(Nt) + 0.5) * (t1 - t0) / Nt
    times_to_emit = np.insert(times_to_emit, 0, 0)  # census particles start at t=0
    angles = -1 + (np.arange(Nu) + 0.5) * 2 / Nu

    NxNu = Nx * Nu  # Number of census particles
    NxNtNu = Nx * Nt * Nu  # Number of source particles

    # Calculate the weight for census particles and source particles
    total_energy_NxNu = a * (T0 ** 4) * Vz
    total_energy_NxNtNu = c * sigma_a * a * (T0 ** 4) * (t1 - t0) * Vz
    weight_census_particle = total_energy_NxNu / NxNu
    weight_source_particle = total_energy_NxNtNu / NxNtNu
    # Store particle information for the deterministic method
    particles = []

    particle_counter = 1  # Counter for particle numbering

    for i in range(Nx):
        for t in times_to_emit:
            for angle in angles:
                particle = {
                    "Name": f"Particle{particle_counter}",
                    "Initial_Position": particle_position[i],
                    "Time_of_Emission": t,
                    "Angle_of_Emission": angle
                }

                # Calculate the distance to census.
                distance_to_census = c * (t1 - particle["Time_of_Emission"]) / abs(angle)

                # Check if census or source particle and then calculate energy deposited
                if t == t0:  # Census particle
                    energy_deposited = weight_census_particle * np.exp(-sigma_a * distance_to_census)

                else:
                    energy_deposited = weight_source_particle * np.exp(-sigma_a * distance_to_census)

                particle["Energy_Deposited"] = energy_deposited
                particle["Distance_to_Census"] = distance_to_census

                particles.append(particle)
                particle_counter += 1

    # Sum energy deposited
    total_energy_deposited = sum(p["Energy_Deposited"] for p in particles)

    return total_energy_deposited


def simulate_monte_carlo(num_particles):
    total_energy_deposited_mc = 0
    seed_value = 5
    np.random.seed(seed_value)

    Nx = num_particles
    Nt = num_particles
    Nu = num_particles
    NxNu = Nx * Nu  # Number of census particles
    NxNtNu = Nx * Nt * Nu  # Number of source particles

    # Calculate the weight for census particles and source particles
    total_energy_NxNu = a * (T0 ** 4) * Vz
    total_energy_NxNtNu = c * sigma_a * a * (T0 ** 4) * (t1 - t0) * Vz
    weight_census_particle = total_energy_NxNu / NxNu
    weight_source_particle = total_energy_NxNtNu / NxNtNu
    # Census particles for the Monte Carlo method
    for i in range(NxNu):
        # Generate random position and angle for the census particle. time = t0
        angle = np.random.uniform(-1, 1)

        # Calculate distance to census and energy deposited for census particle
        distance_to_census_mc = c * (t1 - t0) / abs(angle)  # Distance to census
        energy_deposited = weight_census_particle * np.exp(-sigma_a * distance_to_census_mc)

        # Add the energy deposited to total_energy_deposited
        total_energy_deposited_mc += energy_deposited

    # Source particles for the Monte Carlo method
    for i in range(NxNtNu):
        # Generate random time and angle for the source particle
        time = np.random.uniform(t0, t1)
        angle = np.random.uniform(-1, 1)

        # Calculate distance to census and energy deposited for source particle
        distance_to_census_mc = c * (t1 - time) / abs(angle)  # Distance to census
        energy_deposited = weight_source_particle * np.exp(-sigma_a * distance_to_census_mc)

        # Add the energy deposited to total_energy_deposited
        total_energy_deposited_mc += energy_deposited

    return total_energy_deposited_mc


# (Define Nx, and it assumes Nx=Nt=Nu )
num_particles_to_simulate = np.arange(10, 101, 4, dtype=int)

# # Print Problem Parameters and Results
# print('****************************')
# print('**** Problem Parameters ****')
# print('****************************')
# print()
# print(' Constants used:'
#       '\n Boltzmann constant, k = {} [keV/K]'
#       '\n speed of light, c = {} [cm/ns]'
#       '\n radiation constant, a = {} [keV/cm^3-keV^4]'
#       '\n'
#       '\nOther Properties:'
#       '\nabsorption opacity, sigma_a = {} [cm^-1]'
#       ''.format(k, c, a, sigma_a)
#       )
# print('The zone where particles are sourced is from {} to {} cm.'.format(x0, x1))
# print(f'The initial temperature in the zone is {T0_in_K} K.')
# print('{} particles were generated from {} positions in the source zone,'
#       ' {} time discretizations,'
#       ' and {} discrete angles'
#       ' in a {} ns time step.'.format(Nx * Nt * Nu, Nx, Nt, Nu, t1-t0))
#
# print('Total Energy of initial particles [NxNu] = a * T0^4 * Vz =', total_energy_NxNu)
# print('Total Energy of all source particles [NxNuNt] = c * sigma_a * a * T0^4 * Δt * Vz =', total_energy_NxNtNu)
# print('Ideal Energy deposited = c * sigma_a * a * T0^4 * Δt * Vz=', total_energy_NxNtNu)
# print('Ideal Energy of particles remaining = a * T0^4 * Vz=', total_energy_NxNu)
#
# # Print results for deterministic method
# print()
# print('Results for Deterministic Sourcing Method')
# print()
# print('Energy deposited in the material =', total_energy_deposited)
#
# # Print results for Monte Carlo method
# print()
# print('Results for Monte Carlo Method')
# print()
# print('Energy deposited in the material =', total_energy_deposited_mc)

# Difference from ideal

energy_deposited_deterministic = []
energy_deposited_mc = []

for num_particles in num_particles_to_simulate:
    total_energy_deterministic = simulate_deterministic(num_particles)
    total_energy_mc = simulate_monte_carlo(num_particles)

    energy_deposited_deterministic.append(total_energy_deterministic)
    energy_deposited_mc.append(total_energy_mc)

# Error
ideal_energy_deposited = c * sigma_a * a * (T0 ** 4) * (t1 - t0) * Vz
error_deterministic = np.abs(np.array(energy_deposited_deterministic) - ideal_energy_deposited)
error_mc = np.abs(np.array(energy_deposited_mc) - ideal_energy_deposited)

# Plot
plt.plot(num_particles_to_simulate, error_mc, 'o', label='Monte Carlo')
plt.plot(num_particles_to_simulate, error_deterministic, 'o', label='Deterministic')
plt.legend()
plt.title('Absolute Error of Energy Deposited in 1 Time Step({} ns)'.format(t1))
plt.xlabel('Number of Source Particles (x 10^3)')
plt.ylabel('Absolute Error')
plt.show()
