"""
Title: AEM582 Space Systems Final Project Spring 2025
Authors: Justin Klein, Eli Elstein, Jessica Engel, Misael Alvarez

Created: 03/26/25
Last-Modified: 04/08/25

Inputs: TLE, mass of sc
Outputs:

"""

# imports
import astropy
from astropy import units as u
from astropy.time import Time
import sgp4
from sgp4.api import Satrec
from datetime import datetime
import matplotlib.pyplot as plt
import numpy as np
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.ephem import Ephem
from poliastro.maneuver import Maneuver
from poliastro.plotting import OrbitPlotter
from poliastro.twobody.propagation import CowellPropagator
from poliastro.core.propagation import func_twobody
from poliastro.core.perturbations import J2_perturbation, atmospheric_drag_exponential
from poliastro.constants import rho0_earth, H0_earth
import tletools

starttime = datetime.now()

"""----------------------------------------------------------------------------------"""

# functions
# read in tle
def read_tle_from_file(file_path):
    try:
        with open(file_path, 'r') as file:
            lines = file.readlines()
            if len(lines) < 2:
               raise ValueError("TLE data needs at least two lines.")
            return [line.strip() for line in lines]
    except FileNotFoundError:
        print(f"Error: File not found: {file_path}")
        return None
    except Exception as e:
         print(f"An error occurred: {e}")
         return None
    
# initialize state using TLE
def stateinit(spacecraft):
    # create sgp4 object with tle
    sat = Satrec.twoline2rv(spacecraft.tle1, spacecraft.tle2)

    # find julian date
    jd = sat.jdsatepoch
    fr = sat.jdsatepochF
    startep = Time((jd + fr), format="jd")

    # calculate the spacecraft's position
    e, r, v = sat.sgp4(jd, fr)
    r = r * u.km
    v = v * u.km / u.s
    if e != 0:
        print("Error in TLE")
        return

    # create orbit object
    orb = Orbit.from_vectors(Earth, r, v, epoch=startep)

    # create ephem
    ephem = Ephem.from_orbit(orb, startep)

    # end function
    return ephem, orb

# propogate spacecraft
def propogate(orb, timestep):

    # create propogation time
    epoch = Time(orb.epoch + timestep/86400, format="jd")

    # propogate orbit to next timestep
    neworb = orb.propagate(epoch, method=CowellPropagator(f=j2func))

    # create ephem
    ephem = Ephem.from_orbit(orb, epoch)

    # end function
    return ephem, neworb

# design maneuver
def genburn(spacecraft):
    """ for this I think it would be a good idea to propogate to the next
    apoapsis or periapsis (depending on how the deorbiting occurs) and
    either recircularize at apoapsis or do a hohmann transfer at periapsis
    and then circularize a half an orbit later. So basically only burn the
    most fuel efficient way """
    pass

# J2 pertubations function
def j2func(t0, u_, k):
    du_kep = func_twobody(t0, u_, k)
    ax, ay, az = J2_perturbation(
        t0, u_, k, J2=Earth.J2.value, R=Earth.R.to(u.km).value
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])
    return du_kep + du_ad

# atmospheric drag function
def dragfunc(t0, state, k):
    du_kep = func_twobody(t0, state, k)
    ax, ay, az = atmospheric_drag_exponential(
        t0,
        state,
        k,
        R=R,
        C_D=C_D,
        A_over_m=A_over_m,
        H0=H0,
        rho0=rho0,
    )
    du_ad = np.array([0, 0, 0, ax, ay, az])

    return du_kep + du_ad

"""----------------------------------------------------------------------------------"""

# classes
# spacecraft
class spacecraft:
    def __init__(self, tle1, tle2, mass, fuel, time, burn):
        self.tle1 = tle1
        self.tle2 = tle2
        self.mass = mass
        self.fuel = fuel
        self.time = time
        self.burn = burn
    def print(self):
        print("state:\n", self.tle1, "\n", self.tle2)
        print("mass: ", self.mass)
        print("fuel: ", self.fuel)
        print("time: ", self.time)
        print("burn?: ", self.burn)

"""----------------------------------------------------------------------------------"""

# main code

# initializations
case = 1                            # integer 1-6
tstep = 120                         # seconds
dV = 0                              # current delta V used
dVlimit = 100                       # limit on delta V value
scenario_start_time = 2460774       # julian date (currently April 8, 2025 at 0000z)
mi = 0                              # initial mass
t = list()                          # time since scenario start
statelist = list()                  # list of spacecraft objects at each time step
alt = list()                        # list of spacecraft altitude at each time step
stoploop = False

R = Earth.R.to(u.km).value
k = Earth.k.to(u.km**3 / u.s**2).value
C_D = 2.2                           # coefficient of drag
A_over_m = ((np.pi / 4.0) * (u.m**2) / (100 * u.kg)).to_value(
    u.km**2 / u.kg
)  # km^2/kg
B = C_D * A_over_m
rho0 = rho0_earth.to(u.kg / u.km**3).value  # kg/km^3
H0 = H0_earth.to(u.km).value

# create while loop
while stoploop != True:

    # get spacecraft's current state and ideal state
    if len(t) < 1:
        t.append(0)
        if case == 1:
            tle1 = read_tle_from_file("tles\isstle.txt")[0]
            tle2 = read_tle_from_file("tles\isstle.txt")[1]
            sctle = spacecraft(tle1, tle2, mi, "Hydrogen Peroxide", scenario_start_time, False)
            curephem, curstate = stateinit(sctle)
            idealstate = curstate
            statelist.append(curephem.rv())
            curalt = (np.linalg.norm(curephem.rv()[0]) - R*u.km)/u.km
            alt.append(curalt)
        elif case == 2:
            pass
        elif case == 3:
            pass
        elif case == 4:
            pass
        elif case == 5:
            pass
        elif case == 6:
            pass
    else:
        curephem, curstate = propogate(prevstate, tstep)
        curalt = np.linalg.norm(curephem.rv()[0]) - R*u.km
        if curalt >= 0:
            alt.append(curalt/u.km)
            t.append(t[-1] + tstep)
            statelist.append(curephem.rv())
        else:
            break

    # test if difference of ideal vs current state is above tolerance

    # check if maneuver is currently active

    # if difference is large enough and maneuver is currently not active, design maneuver

    # calculate Delta V usage

    # determine if Delta V limit has been reached

    # amend time list (and any other performance lists)

    # move to next time step or exit while loop
    if dV >= dVlimit:
        stoploop = True
    else:
        prevstate = curstate
    
    if len(t) == 72000:
        stoploop = True

# plot outputs
plt.figure()
plt.plot(t, alt)

# export output data if necessary

# exit main code

"""----------------------------------------------------------------------------------"""

# output run time and show plots
endtime = datetime.now()
print("Runtime is", endtime-starttime)
plt.show()