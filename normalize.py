import sys
import csv
import math

scans = {} # empty dictionary.  Will contain one entry per scan.
telomere_names=[]
with open(sys.argv[1]) as csvfile:
    b_reader = csv.reader(csvfile)
    row_count = 0
    for row in b_reader:
        if row_count == 0:
            row_pos = 0
            for str in row:
                if row_pos>0 and (row_pos-1) % 5 == 0:
                    telomere_names.append(str.split("_")[0])
                row_pos += 1
        else:
            scan_values = {}
            scans[row[0]] = scan_values # will contain a dictionary of telomeres
            row_pos = 0
            for str in row:
                t_name = ""
                if row_pos>0 and (row_pos-1) % 5 == 0:
                    t_name = telomere_names[int((row_pos-1)/5)]
                    t_values = {}
                    scan_values[t_name] = t_values # will contain all descriptors of this telomere
                    t_values["x"] = float(str)
                elif row_pos>1 and (row_pos-2) % 5 == 0:
                    t_values["y"] = float(str)
                row_pos += 1
        row_count += 1
#  print(scans)

for scankey in scans:
    scan = scans[scankey]
    for str in ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]:
        scan.pop(str)

    for t_mere_key in scan:
        t_mere = scan[t_mere_key]

        # t_mere["z1"] = math.sqrt(25 - t_mere["x"]**2 - t_mere["y"]**2)
        # t_mere["z2"] = -t_mere["z1"]
        if t_mere["x"]**2 + t_mere["y"]**2 > 25:
            print(scankey)
            print(t_mere_key)
            print(t_mere)



