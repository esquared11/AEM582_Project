"""
Title: AEM582 Space Systems Final Project Spring 2025
Authors: Justin Klein, Eli Elstein, Jessica Engel, Misael Alvarez

Created: 03/26/25
Last-Modified: 04/18/25

Inputs: TLE
Outputs: plots displaying life-of-mission altitude, etc.

"""

# imports
import os
from datetime import datetime
from sgp4.api import Satrec
from hapsira.core.elements import eccentricity_vector
from hapsira.bodies import Earth
from hapsira.twobody import Orbit
from hapsira.ephem import Ephem
from astropy import units as u
from astropy.time import Time
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks
from scipy.signal import argrelextrema

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
def stateinit(tle):
    # create sgp4 object with tle
    sat = Satrec.twoline2rv(tle[0], tle[1])

    # find julian date
    jd = sat.jdsatepoch
    fr = sat.jdsatepochF
    startep = Time((jd + fr), format="jd")

    # calculate the spacecraft's position
    e, r, v = sat.sgp4(jd, fr)

    # calculate eccentricity
    k = Earth.k.to(u.km**3 / u.s**2).value
    eccv = eccentricity_vector(k, np.array(r), np.array(v))
    ecc = np.linalg.norm(eccv)

    # calculate periapsis and apoapsis
    orb = Orbit.from_vectors(Earth, (r*u.km), (v*u.km/u.s), epoch=startep)

    return [r, v], ecc, orb.r_p/u.km, orb.r_a/u.km

# main code
starttime = datetime.now()

# read in tle
tle = read_tle_from_file("tles\\isstle.txt")
state, e, r_p, r_a = stateinit(tle)
if e < 0.001:
    desiredrmag = int(np.linalg.norm(state[0]))
    altthresh = int(desiredrmag - 6405)
else:
    desiredrmag = int(r_a)
    altthresh = int(r_a - 6410)
    deorbitthresh = altthresh - 180

# append gmat script file
with open("test.script", 'r') as gmatfile:
    data = gmatfile.readlines()

data[11] = "LEOsat.X = " + str(state[0][0]) + "\n"
data[12] = "LEOsat.Y = " + str(state[0][1]) + "\n"
data[13] = "LEOsat.Z = " + str(state[0][2]) + "\n"
data[14] = "LEOsat.VX = " + str(state[1][0]) + "\n"
data[15] = "LEOsat.VY = " + str(state[1][1]) + "\n"
data[16] = "LEOsat.VZ = " + str(state[1][2]) + "\n"
data[71] = "desiredRMAG = " + str(desiredrmag) + "\n"
data[72] = "desiredECC = " + str(e) + "\n"

if e < 0.001:
    data[87] = "Propagate 'Prop One Step' LEOprop(LEOsat) \n"
    data[88] = "\n"
    data[90] = "If 'If Alt < Threshold' LEOsat.Earth.Altitude < " + str(altthresh) + "\n"
    data[113] = "\n"
    data[114] = "\n"
    data[115] = "Propagate LEOprop(LEOsat) {LEOsat.Earth.Altitude = 250}\n"
    data[116] = "\n"
    data[117] = "\n"
else:
    data[87] = "Propagate 'Prop To Apoapsis' LEOprop(LEOsat) {LEOsat.Apoapsis} \n"
    data[88] = "\n"
    data[90] = "If 'If Alt < Threshold' LEOsat.Earth.Altitude < " + str(altthresh) + "\n"
    data[113] = "scapoalt = 50000000 \n"
    data[114] = "While scapoalt > " + str(deorbitthresh) + "\n"
    data[115] = "Propagate LEOprop(LEOsat) {LEOsat.Apoapsis} \n"
    data[116] = "scapoalt = LEOsat.Earth.Altitude \n"
    data[117] = "EndWhile \n"

with open("test.script", 'w') as gmatfile:
    gmatfile.writelines(data)

gmatfile.close()

# create system call and run
args = [r"C:\\Users\\eelstein\\GMAT\\bin\\GMAT.exe","--logfile",
        r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\gmatlog.txt",
        "--run", r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\test.script"]

test = r' '.join(args)

os.system(test)


# open output files
filename = "C:\\Users\\eelstein\\GMAT\\output\\rf2.txt"
file = open(filename)

# get data from output files
timelist = list()
altlist = list()
rmaglist = list()
ecclist = list()
deltav = float()
inclist = list()
curdv = 0
prevdv = 0
for i, line in enumerate(file):
    if i == 0:
        pass
    elif i == 1:
        newline = line.split()
        startday = float(newline[0])
        prevdv = abs(float(newline[4])) + abs(float(newline[5]))
    else:
        newline = line.split()
        timelist.append(float(newline[0]) - startday)
        altlist.append(float(newline[1]))
        rmaglist.append(float(newline[2]))
        ecclist.append(float(newline[3]))
        curdv = abs(float(newline[4])) + abs(float(newline[5])) + abs(float(newline[8]))
        if curdv != prevdv:
            deltav += curdv
            prevdv = curdv
        inclist.append(float(newline[7]))

# thin altitude data into periapsis and apoapsis
altnp = np.array(altlist)
timenp = np.array(timelist)
peaks, _ = find_peaks(altlist)
mins = argrelextrema(altnp, np.less)

# plot data
plt.figure(1)
plt.plot(timelist, altlist, 'dimgrey', label='Altitude')
plt.plot(timenp[peaks], altnp[peaks], 'r.', label='Apoapsis')
plt.plot(timenp[mins], altnp[mins], 'b.', label='Periapsis')
plt.xlabel("Time (Days)")
plt.ylabel("Altitude (km)")
plt.title("Case 1 Altitude")
plt.legend()
plt.figure(2)
plt.plot(timelist, ecclist, 'dimgrey')
plt.xlabel("Time (Days)")
plt.ylabel("Eccentricity")
plt.title("Case 1 Eccentricity")
plt.figure(3)
plt.plot(timelist, inclist, 'dimgrey')
plt.xlabel("Time (Days)")
plt.ylabel("Inclination (Deg)")
plt.title("Case 1 Inclination")

totime = timelist[-1] - timelist[0]

print("Total Delta V:", deltav*1000, "m/s")
print("Total Time:", totime, "(Days)")

file.close()

# end script
endtime = datetime.now()
print("Run Time: ", endtime-starttime)

plt.show()