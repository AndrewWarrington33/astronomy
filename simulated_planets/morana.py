# Modeling the theoretical solar system Morana with the Rebound package
import numpy as np
import matplotlib.pyplot as plt
import rebound
from astropy import units as u
import time

sim = rebound.Simulation()
sim.units = ('yr', 'AU', 'Msun')
sim.add(m = 0.867)
sim.add(m = (1*u.Mearth).to(u.Msun).value,
        a = 0.7,
        e =0.5,
        omega = 6.01)
sim.add(m = (5*u.Mearth).to(u.Msun).value,
        a = 2.1,
        e = 0.1,
        omega = 2.85,
        inc = (6.5*u.deg).to(u.rad).value,
        Omega = 4.05)
sim.add(m = (93*u.Mearth).to(u.Msun).value,
        a = 6.3,
        e = 0.3,
        omega = 5.86,
        inc = (21*u.deg).to(u.rad).value,
        Omega = 5.72)
sim.add(m=(.6*u.Mearth).to(u.Msun).value,
        a = .51,
        e = 0.12,
        omega = 5.37,
        inc = (-12*u.deg).to(u.rad).value,
        Omega = 2.14)
sim.add(primary=sim.particles[1],
        m = (2.1e17*u.kg).to(u.Msun).value,
        a = (67000*u.km).to(u.AU).value,
        e = 0.09,
        omega = 2.35,
        inc = (5.2*u.deg).to(u.rad).value,
        Omega = 4.87)
sim.add(primary=sim.particles[1],
        m = (4.7e15*u.kg).to(u.Msun).value,
        a = (150000*u.km).to(u.AU).value,
        e = 0.23,
        omega = 3.78,
        inc = (-17*u.deg).to(u.rad).value,
        Omega = 3.98)
sim.add(primary=sim.particles[1],
        m = (6.1e20*u.kg).to(u.Msun).value,
        a = (232000*u.km).to(u.AU).value,
        e = 0.04,
        omega = 5.10,
        inc = (1.4*u.deg).to(u.rad).value,
        Omega = 0.77)

fig = rebound.OrbitPlot(sim)
op1 = rebound.OrbitPlot(sim, color=True, particles = [1,2,3,4,])
op2 = rebound.OrbitPlot(sim, color=True, particles = [5,6,7], primary=1, show_primary=True)
plt.show()

p = sim.particles

times = np.linspace(0, 1000)
smas = np.full((len(p)-1,len(times)), np.nan)
ds = np.zeros((len(p)-1,len(times)))
time0 = time.time()
for i, t in enumerate(times):
        print(t)
        star_loc = [p[0].x, p[0].y, p[0].z]
        for j in range(1,len(p)):
                ds[j-1,i] = np.sqrt((p[j].x-star_loc[0])**2 
                                    + (p[j].y-star_loc[1])**2 + 
                                    (p[j].z-star_loc[2])**2)
                if j <=4:
                        smas[j-1,i] = p[j].a
                else:
                        orb = p[j].orbit(primary=p[1])     
                        smas[j-1,i] = orb.a 
        sim.integrate(t)
        
print(time.time()-time0, flush=True)

plt.figure()
plt.plot(times,smas.T)

plt.figure()
plt.plot(times,ds.T)
plt.show()