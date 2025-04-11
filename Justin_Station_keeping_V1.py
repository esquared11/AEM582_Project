import numpy as np
import matplotlib.pyplot as plt
from tqdm import tqdm
from datetime import timedelta

from sgp4.api import Satrec, jday
from astropy import units as u
from astropy.time import Time, TimeDelta
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

# ---- Configuration Parameters ----
DECAY_THRESHOLD = 0.02                # Deviation threshold for station-keeping
ISP = 300                             # Specific impulse in seconds
G0 = 9.80665                          # Standard gravity [m/s²]
INITIAL_PROPELLANT_KG = 50.0         # Starting propellant
DRY_MASS_KG = 500.0                  # Dry mass of satellite
REENTRY_ALTITUDE = 120 * u.km        # Reentry threshold
BURN_INTERVAL = 12 * u.hour          # Min time between burns
STEP = 1 * u.hour                    # Propagation step size
MAX_DURATION = 5 * 365 * u.day       # Optional simulation duration cap

# ---- Hardcoded TLEs ----
TLES = [
    {
        "name": "ISS",
        "tle_line1": "1 25544U 98067A   24086.47376251  .00012171  00000+0  22009-3 0  9993",
        "tle_line2": "2 25544  51.6409  24.0604 0001747 110.9241 313.4433 15.50750418393593"
    },
    {
        "name": "NOAA 19",
        "tle_line1": "1 33591U 09005A   24086.55222955  .00000099  00000+0  86047-4 0  9997",
        "tle_line2": "2 33591  99.1956  25.7748 0014191  80.1248 280.2186 14.12420987799964"
    }
]

# ---- Main Simulation Function ----
def simulate_satellite(name, line1, line2):
    print(f"\n Starting simulation for {name}...")

    # Initialize orbit using SGP4
    sat = Satrec.twoline2rv(line1, line2)
    tle_epoch = sat.jdsatepoch + sat.jdsatepochF
    dt_start = Time(tle_epoch, format="jd")

    jd, fr = jday(dt_start.datetime.year, dt_start.datetime.month, dt_start.datetime.day,
                  dt_start.datetime.hour, dt_start.datetime.minute, dt_start.datetime.second)
    e, r, v = sat.sgp4(jd, fr)
    if e != 0:
        print(f"Error propagating {name} from TLE.")
        return

    r = np.array(r) * u.km
    v = np.array(v) * u.km / u.s
    orbit = Orbit.from_vectors(Earth, r, v, epoch=dt_start)

    # Initial tracking state
    propellant_remaining = INITIAL_PROPELLANT_KG
    last_burn_time = orbit.epoch - BURN_INTERVAL
    initial_altitude = np.linalg.norm(orbit.r) - Earth.R
    initial_inclination = orbit.inc
    fuel_depleted = False
    decayed = False

    # Logging arrays
    epochs, altitudes, inclinations, dv_log, fuel_log = [], [], [], [], []

    current_orbit = orbit
    steps_estimate = int((MAX_DURATION / STEP).to_value(u.dimensionless_unscaled))

    with tqdm(total=steps_estimate, desc=f"{name}", unit="step") as pbar:
        while not decayed:
            alt = np.linalg.norm(current_orbit.r) - Earth.R
            inc = current_orbit.inc

            # Log data
            epochs.append(current_orbit.epoch.datetime)
            altitudes.append(alt.to_value(u.km))
            inclinations.append(inc.to_value(u.deg))
            fuel_log.append(propellant_remaining)

            # Deviation checks
            alt_dev = abs(alt - initial_altitude) / initial_altitude
            inc_dev = abs(inc - initial_inclination) / initial_inclination
            out_of_bounds = (alt_dev > DECAY_THRESHOLD) or (inc_dev > DECAY_THRESHOLD)
            time_since_burn = (current_orbit.epoch - last_burn_time).to(u.hour)

            # Burn logic
            if not fuel_depleted and out_of_bounds and time_since_burn >= BURN_INTERVAL:
                orb1 = current_orbit
                orb2 = Orbit.circular(Earth, alt=initial_altitude, inc=initial_inclination)
                man = Maneuver.hohmann(orb1, orb2.a)
                dv = man.get_total_cost().to_value(u.m / u.s)

                m0 = DRY_MASS_KG + propellant_remaining
                mf = m0 / np.exp(dv / (ISP * G0))
                fuel_used = m0 - mf

                if fuel_used <= propellant_remaining:
                    propellant_remaining -= fuel_used
                    last_burn_time = current_orbit.epoch
                    dv_log.append(dv)
                else:
                    fuel_depleted = True
                    dv_log.append(0)
            else:
                dv_log.append(0)

            # Reentry or time cap
            if alt < REENTRY_ALTITUDE:
                print(f"{name} decayed at {current_orbit.epoch.datetime}")
                break

            if (current_orbit.epoch - dt_start) > MAX_DURATION:
                print(f"{name} hit simulation cap of {MAX_DURATION.to(u.day):.0f}")
                break

            # Advance orbit
            current_orbit = current_orbit.propagate(STEP)
            pbar.update(1)

    # Final report
    mission_duration = current_orbit.epoch - dt_start
    print(f"{name} lifetime: {mission_duration.to(u.day):.2f}")

    # Plot results
    plt.figure(figsize=(14, 9))
    plt.suptitle(f"{name} Lifetime Propagation", fontsize=14)

    plt.subplot(4, 1, 1)
    plt.plot(epochs, altitudes)
    plt.axhline(initial_altitude.to_value(u.km) * (1 - DECAY_THRESHOLD), color='r', linestyle='--')
    plt.axhline(initial_altitude.to_value(u.km) * (1 + DECAY_THRESHOLD), color='r', linestyle='--')
    plt.ylabel("Altitude [km]")
    plt.grid(True)

    plt.subplot(4, 1, 2)
    plt.plot(epochs, inclinations, color='purple')
    plt.axhline(initial_inclination.to_value(u.deg) * (1 - DECAY_THRESHOLD), color='r', linestyle='--')
    plt.axhline(initial_inclination.to_value(u.deg) * (1 + DECAY_THRESHOLD), color='r', linestyle='--')
    plt.ylabel("Inclination [deg]")
    plt.grid(True)

    plt.subplot(4, 1, 3)
    plt.plot(epochs, fuel_log, color='green')
    plt.ylabel("Fuel [kg]")
    plt.grid(True)

    plt.subplot(4, 1, 4)
    plt.plot(epochs, dv_log, color='orange')
    plt.xlabel("UTC Time")
    plt.ylabel("ΔV [m/s]")
    plt.grid(True)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    plt.show()

# ---- Run All Simulations ----
if __name__ == "__main__":
    for sat in TLES:
        simulate_satellite(sat["name"], sat["tle_line1"], sat["tle_line2"])
