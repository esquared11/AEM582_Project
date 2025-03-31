"""
Title: AEM582 Space Systems Final Project Spring 2025
Authors: Jessica Engel, Misael Alvarez, Justin Klein, Eli Elstein

Created: 03/26/25
Last-Modified: 03/26/25

Inputs:
Outputs:

"""

# imports
import astropy
import sgp4
from datetime import datetime
#import matplotlib.pyplot as plt

import numpy
from astropy import units as u
from poliastro.bodies import Earth
from poliastro.twobody import Orbit
from poliastro.maneuver import Maneuver

starttime = datetime.now()

# main code

#TLE Read in - Justin
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






# output run time
endtime = datetime.now()
print("Runtime is", endtime-starttime)
