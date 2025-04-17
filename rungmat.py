"""

"""

# imports
import os
from datetime import datetime

# main code
starttime = datetime.now()

args = [r"C:\\Users\\eelstein\\GMAT\\bin\\GMAT.exe",
        "--logfile", r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\gmatlog.txt",
        "--run", r"C:\\Users\\eelstein\\Documents\\Bama\\AEM582\\AEM582_Project\\test.script"]

test = r' '.join(args)

os.system(test)

endtime = datetime.now()

print("Run Time: ", endtime-starttime)