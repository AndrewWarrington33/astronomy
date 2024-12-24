# Modeling the system of bodies with the Rebound package
import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
import time

# Create the simulation
sim = rebound.Simulation()

# Adding bodies to the simulation (non-solar system objects)
sim.add(m = 0.957)  # Primary body
sim.add(m = 0.342,
        a = 0.08145,
        e = 0.0288,
        inc = np.deg2rad(90-89.613),
        omega = np.deg2rad((226.3)))  # First orbiting body

sim.add(m = (2.07*u.Mearth).to(u.Msun).value,
        a = .2877,
        e = .021,
        inc = np.deg2rad(90-89.752),
        omega = np.deg2rad(48.6))  # Second orbiting body

sim.add(m = (19.02*u.Mearth).to(u.Msun).value,
        a = .6992,
        e = .041,
        inc = np.deg2rad(90-90.395),
        omega = np.deg2rad(352))  # Third orbiting body

sim.add(m = (3.17*u.Mearth).to(u.Msun).value,
        a = .9638,
        e = .044,
        inc = np.deg2rad(90-90.1925),
        omega = np.deg2rad(306))  # Fourth orbiting body

# Adding a moon-like body (satellite of particle 4)
sim.add(primary = sim.particles[4],
         P = (5.877*u.day).to(u.yr).value,
         e = 0.01,
         omega = 4.77,
         inc = np.deg2rad(4),
         Omega = 0.83)

# Moving the system to the center of mass
sim.move_to_com()

# Plotting orbits of bodies in the system
op1 = rebound.OrbitPlot(sim, color=True, particles=[1, 2, 3, 4])
op2 = rebound.OrbitPlot(sim, color=True, particles=[5], primary=4, show_primary=True)
plt.show()

# Simulating the orbits over time
p = sim.particles
times = np.linspace(0, 12/365, 1000)  # Time period in years
smas = np.full((len(p)-1, len(times)), np.nan)  # Semi-major axes array
ds = np.full((len(p)-1, len(times)), np.nan)  # Distance array

time0 = time.time()

# Loop to integrate the simulation over the time steps
for i, t in enumerate(times):
    print(f"Integrating time step {i+1}/{len(times)}: {t} years")
    
    # Save the star location (primary body)
    star_loc = [p[0].x, p[0].y, p[0].z]
    plt.plot(p[0].x, p[0].y, '.', color='C0', alpha=i/len(times))  # Primary body
    
    # Plot the positions of orbiting bodies
    for j in range(1, len(p)):
        plt.plot(p[j].x, p[j].y, '.', color='C'+str(j), alpha=i/len(times))
        ds[j-1, i] = np.sqrt((p[j].x - star_loc[0])**2 + 
                             (p[j].y - star_loc[1])**2 + 
                             (p[j].z - star_loc[2])**2)  # Distance from primary body
        
        # Update semi-major axes for bodies
        if j <= 4:
            smas[j-1, i] = p[j].a  # Semi-major axis for primary orbiting bodies
        else:
            orb = p[j].orbit(primary=p[1])  # Orbit of moon-like body around particle 1
            smas[j-1, i] = orb.a

    # Integrate the simulation by one time step
    sim.integrate(t)

# Set plot limits
lim = 1.1
plt.xlim(-lim, lim)
plt.ylim(-lim, lim)
plt.show()

# Print simulation time
print(f"Simulation time: {time.time()-time0} seconds", flush=True)

# Plot semi-major axes over time
plt.figure()
plt.plot(times, smas.T)
plt.title('Semi-Major Axes Over Time')
plt.xlabel('Time (years)')
plt.ylabel('Semi-Major Axis (AU)')

# Plot distances over time
plt.figure()
plt.plot(times, ds.T)
plt.title('Distances from Primary Body Over Time')
plt.xlabel('Time (years)')
plt.ylabel('Distance (AU)')

plt.show()
