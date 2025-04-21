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
import numpy as np
from matplotlib import pyplot as plt
from scipy.signal import find_peaks

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

    # calculate the spacecraft's position
    e, r, v = sat.sgp4(jd, fr)

    return [r, v]

# main code
starttime = datetime.now()

# read in tle
tle = read_tle_from_file("tles\isstle.txt")
state = stateinit(tle)
desiredrmag = int(np.linalg.norm(state[0]))
altthresh = int(desiredrmag - 6385)

# append gmat script file
with open("test.script", 'r') as gmatfile:
    data = gmatfile.readlines()

data[11] = "LEOsat.X = " + str(state[0][0]) + "\n"
data[12] = "LEOsat.Y = " + str(state[0][1]) + "\n"
data[13] = "LEOsat.Z = " + str(state[0][2]) + "\n"
data[14] = "LEOsat.VX = " + str(state[1][0]) + "\n"
data[15] = "LEOsat.VY = " + str(state[1][1]) + "\n"
data[16] = "LEOsat.VZ = " + str(state[1][2]) + "\n"
data[68] = "desiredRMAG = " + str(desiredrmag) + "\n"
data[83] = "If 'If Alt < Threshold' LEOsat.Earth.Altitude < " + str(altthresh) + "\n"

with open("test.script", 'w') as gmatfile:
    gmatfile.writelines(data)

gmatfile.close()

# create system call and run
args = [r"C:\\Users\\eelstein\\GMAT\\bin\\GMAT.exe",
        "--logfile", r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\gmatlog.txt",
        "--run", r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\test.script"]

test = r' '.join(args)

os.system(test)

# wait for gmat completion and then exit out of gmat





# open output files
filename = "C:\\Users\\eelstein\\GMAT\\output\\rf2.txt"
file = open(filename)

# get data from output files
timelist = list()
altlist = list()
rmaglist = list()
ecclist = list()
deltav = float()
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
        curdv = abs(float(newline[4])) + abs(float(newline[5]))
        if curdv != prevdv:
            deltav += curdv
            prevdv = curdv

# thin altitude data into periapsis and apoapsis
altnp = np.array(altlist)
timenp = np.array(timelist)
peaks, _ = find_peaks(altlist)

# plot data
plt.figure(1)
plt.plot(timelist, altlist)
plt.plot(timenp[peaks], altnp[peaks], 'ro')
plt.xlabel("Time (Days)")
plt.ylabel("Altitude (km)")
plt.figure(2)
plt.plot(timelist, ecclist)
plt.xlabel("Time (Days)")
plt.ylabel("Eccentricity")
plt.show()

print("Total Delta V:", deltav*1000, "m/s")

file.close()

# end script
endtime = datetime.now()
print("Run Time: ", endtime-starttime)