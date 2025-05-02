# calculate_stationkeeping_dv.py

# Imports
import numpy as np
import matplotlib.pyplot as plt

# --- Step 1: Read the GMAT output file ---
def read_gmat_output(filepath):
    timelist = []
    altlist = []
    rmaglist = []
    ecclist = []
    inclist = []
    deltavlist = []
    deltav = 0
    deltav_cumulative = 0.0
    curdv = 0.0
    prevdv = 0.0
    burnnum = 0
    instant_deltav = []
    velx = list()
    vely = list()
    velz = list()

    with open(filepath, 'r') as file:
        for i, line in enumerate(file):
            if i == 0:
                continue
            elif i == 1:
                newline = line.split()
                startday = float(newline[0])
                prevdv = abs(float(newline[4])) + abs(float(newline[5]))
            else:
                newline = line.split()
                time = float(newline[0]) - startday
                alt = float(newline[1])
                rmag = float(newline[2])
                ecc = float(newline[3])
                curdv = abs(float(newline[4])) + abs(float(newline[5])) + abs(float(newline[8]))
                deltav_change = 0.0
                if curdv != prevdv:
                    deltav += curdv
                    deltav_change = (curdv - prevdv)
                    instdeltv = curdv
                    deltav_cumulative += deltav_change
                    prevdv = curdv
                    burnnum += 1
                else:
                    instdeltv = 0
                incl = float(newline[7])

                timelist.append(time)
                altlist.append(alt)
                rmaglist.append(rmag)
                ecclist.append(ecc)
                inclist.append(incl)
                deltavlist.append(deltav)
                instant_deltav.append(instdeltv)

                velx.append(float(newline[9]))
                vely.append(float(newline[10]))
                velz.append(float(newline[11]))

    return (np.array(timelist), np.array(altlist), np.array(rmaglist), 
            np.array(ecclist), np.array(inclist), 
            np.array(deltavlist), np.array(instant_deltav), velx, vely, velz, burnnum)


# --- Step 2: Set the correct output file path ---
output_file = "C:\\Users\\eelstein\\GMAT\\output\\rf2.txt"  # Adjust path if needed, I kept it to you system's path for this update

# --- Step 3: Run the function ---
(timelist, altlist, rmaglist, ecclist, inclist, 
 deltavlist, instant_deltavlist, velx, vely, velz, burnnum) = read_gmat_output(output_file)

# --- Step 4: Report the Results ---
print(f"Total Station-Keeping Delta-V Required: {deltavlist[-1] * 1000:.2f} m/s")
print(f"Total Time Simulated: {timelist[-1]:.2f} days")
avburnstr = deltavlist[-1]/burnnum
print("Average Burn Strength: ", avburnstr)

# --- Step 5: Plot Altitude vs. Time ---
plt.figure()
plt.plot(timelist, altlist, color='black')
plt.title("Altitude vs. Time")
plt.xlabel("Time (days)")
plt.ylabel("Altitude (km)")
plt.grid(True)

# --- Step 6: Plot Cumulative Delta-V vs. Time ---
plt.figure()
plt.plot(timelist, deltavlist * 1000, color='blue')
plt.title("Cumulative Delta-V vs. Time")
plt.xlabel("Time (days)")
plt.ylabel("Cumulative Delta-V (m/s)")
plt.grid(True)

# --- Step 7: Plot Instantaneous Delta-V Usage ---
plt.figure()
plt.plot(timelist, instant_deltavlist * 1000, color='red')
plt.title("Instantaneous Delta-V Usage vs. Time")
plt.xlabel("Time (days)")
plt.ylabel("Instantaneous Delta-V (m/s)")
plt.grid(True)


plt.figure()
plt.plot(timelist, velx, '--', label='VX')
plt.plot(timelist, vely, '--', label='VY')
plt.plot(timelist, velz, '--', label='VZ')
plt.title("Spacecraft Veloicty by Components")
plt.xlabel("Time (days)")
plt.ylabel("Velocity (km/s)")
plt.legend()

plt.show()
