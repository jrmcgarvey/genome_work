import math
import random
import csv
import sys
import utils


outfile=sys.argv[1]
scancount = int(sys.argv[2])
t_names = ["Y1p", "Xp", "Xq", "Y1q", "Y2p", "Y2q"]
scans={}


scan={}
for t_name in t_names:
    t_mere = {}
    scan[t_name] = t_mere 
    t_mere["p"]=utils.random_point()
    
scans["scan0"] = scan

for i in range(1,scancount):
    r = utils.random_rotation()
    scani={}
    scans[f"scan{i}"] = scani
    for t_name,t_mere in scan.items():
        new_t_mere={}
        scani[t_name] = new_t_mere
        new_t_mere["p"] = utils.matrix_rotate(r,t_mere["p"])

with open(sys.argv[1],"w") as csvfile:
    twriter = csv.writer(csvfile)
    row = [""]
    for t_name in t_names:
        row.extend( [t_name + "_x",t_name + "_y","z","",""])
    twriter.writerow(row)
    for scankey, scan in scans.items():
        row = [scankey]
        for t_name, t_mere in scan.items():
            row.extend([t_mere["p"][0], t_mere["p"][1],t_mere["p"][2],"",""])
        twriter.writerow(row)

