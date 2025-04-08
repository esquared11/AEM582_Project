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
import sgp4
from sgp4.api import Satrec
from datetime import datetime
#import matplotlib.pyplot as plt
import numpy
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

starttime = datetime.now()

"""----------------------------------------------------------------------------------"""

# functions
# read in TLE - Justin
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

# propogate spacecraft
def propogate(spacecraft, timestep):
    # check if burn is detected
    if spacecraft.burn == False:
        # create sgp4 object with spacecraft state
        sat = Satrec.twoline2rv(spacecraft.tle1, spacecraft.tle2)

        # find julian date
        jd = int(spacecraft.time + timestep) 
        fr = spacecraft.time + timestep - jd

        # calculate the spacecraft's position
        e, r, v = sat.sgp4(jd, fr)

        # generate new tle
        #new_spacecraft = newtle(e, r, v)        dont actually need a function just do it here

    elif spacecraft.burn == True:
        pass # put propogation with maneuver in here

    # end function and return new state
    #return new_spacecraft

# design maneuver
def genburn(spacecraft):
    """ for this I think it would be a good idea to propogate to the next
    apoapsis or periapsis (depending on how the deorbiting occurs) and
    either recircularize at apoapsis or do a hohmann transfer at periapsis
    and then circularize a half an orbit later. So basically only burn the
    most fuel efficient way """
    pass

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
        print("state:\n", self.tle1, "\n ", self.tle2)
        print("mass: ", self.mass)
        print("fuel: ", self.fuel)
        print("time: ", self.time)
        print("burn?: ", self.burn)

"""----------------------------------------------------------------------------------"""

# main code

# initializations
case = 1                            # integer 1-6
tstep = 1                           # seconds
dV = 0                              # current delta V used
dVlimit = 100                       # limit on delta V value
scenario_start_time = 2460774       # julian date (currently April 8, 2025 at 0000z)
t = list()                          # time since scenario start
statelist = list()                  # list of spacecraft objects at each time step
stoploop = False

# create while loop
while stoploop != True:

    # get spacecraft's current state
    if len(t) < 1:
        t.append(0)
        if case == 1:
            pass
            # this will be grabbing a specific TLE and fuel stats, something like curstate = blah
            curstate = 'something'
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
        t.append(t[-1] + 1)
        # propogate from latest state
        curstate = propogate(prevstate, tstep)

    # get spacecraft's ideal state

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
        curstate = 'nothing'
    #if len(t) == 9000000:
    #    stoploop = True

# exit main code

"""----------------------------------------------------------------------------------"""

# output run time
endtime = datetime.now()
print("Runtime is", endtime-starttime)
