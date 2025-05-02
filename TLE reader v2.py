from datetime import datetime, timedelta
import numpy as np

def parse_tle(tle_line1, tle_line2):
    # Extract values from fixed-width TLE fields
    epoch_year = int(tle_line1[18:20])
    epoch_day = float(tle_line1[20:32])
    mean_motion = float(tle_line2[52:63])
    
    inclination = float(tle_line2[8:16])
    raan = float(tle_line2[17:25])
    eccentricity = float('0.' + tle_line2[26:33])
    arg_perigee = float(tle_line2[34:42])
    mean_anomaly = float(tle_line2[43:51])
    
    # Convert epoch to datetime
    year = 2000 + epoch_year if epoch_year < 57 else 1900 + epoch_year
    epoch = datetime(year, 1, 1) + timedelta(days=epoch_day - 1)

    # Compute semi-major axis (a) from mean motion (n)
    mu = 398600.4418  # km^3/s^2 (Earth)
    n_rad = mean_motion * 2 * np.pi / 86400  # rev/day -> rad/s
    a = (mu / n_rad**2)**(1/3)

    # Compute perigee/apogee
    perigee = a * (1 - eccentricity) - 6371  # altitude above Earth surface
    apogee = a * (1 + eccentricity) - 6371
    period = 86400 / mean_motion  # seconds

    return {
        "Epoch (UTC)": epoch,
        "Inclination (deg)": round(inclination, 4),
        "RAAN (deg)": round(raan, 4),
        "Eccentricity": round(eccentricity, 7),
        "Argument of Perigee (deg)": round(arg_perigee, 4),
        "Mean Anomaly (deg)": round(mean_anomaly, 4),
        "Mean Motion (rev/day)": round(mean_motion, 8),
        "Semi-major Axis (km)": round(a, 2),
        "Perigee Altitude (km)": round(perigee, 2),
        "Apogee Altitude (km)": round(apogee, 2),
        "Orbital Period (min)": round(period / 60, 2),
    }

# ISS TLE
tle_line1 = "1 25544U 98067A   24086.47376251  .00012171  00000+0  22009-3 0  9993"
tle_line2 = "2 25544  51.6409  24.0604 0001747 110.9241 313.4433 15.50750418393593"

# JD 3 TLE
# tle_line1 = "1 14795U 84012F   25118.76645763  .00010014  00000-0  74922-3 0  9995"
# tle_line2 = "2 14795  63.3010  40.8942 0770515 118.9436 249.1100 13.83738910297182"

# Run parser
tle_data = parse_tle(tle_line1, tle_line2)

# Display results
for key, value in tle_data.items():
    print(f"{key}: {value}")