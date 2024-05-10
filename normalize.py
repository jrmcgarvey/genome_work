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

short_distances = {}
long_distances = {}

for t_name in ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]:
    telomere_names.remove(t_name)

for t_name in telomere_names:
    if t_name != "Y1p":
        short_distances[t_name]=[]
        long_distances[t_name]=[]  # these arrays will help us identify plausible distances

for scankey in scans:
    scan = scans[scankey]
    for str in ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]:
        if scan.get(str):
            scan.pop(str)
    bad_values = []
    for t_mere_key in scan:
        t_mere = scan[t_mere_key]

        if t_mere["x"]**2 + t_mere["y"]**2 > 25:
            bad_values.append(t_mere_key) # throw out the bad values for now
        else:
            t_mere["z1"] = math.sqrt(25 - t_mere["x"]**2 - t_mere["y"]**2)
            t_mere["z2"] = -t_mere["z1"]
    for bv in bad_values:
        scan.pop(bv)
    ref_point = scan["Y1p"]
    ref_point.pop("z2") # assume z value for Y1p is positive
    for t_mere_key in scan:
        if t_mere_key != "Y1p":
            t_mere = scan[t_mere_key]
            t_mere["dist1"] = (ref_point["x"] - t_mere["x"])**2 + (ref_point["y"] - t_mere["y"])**2 + \
                              (ref_point["z1"] - t_mere["z1"])**2
            t_mere["dist2"] =  (ref_point["x"] - t_mere["x"])**2 + (ref_point["y"] - t_mere["y"])**2 + \
                              (ref_point["z1"] - t_mere["z2"])**2
            short_distances[t_mere_key].append((scankey,t_mere["dist1"]))
            long_distances[t_mere_key].append((scankey,t_mere["dist2"]))

# Ok, now we've computed all the distances, let's throw away the ones that aren't probable
# print(short_distances)
# print(long_distances)

def key_for_sort(e):
    return e[1]
    
for t_name in telomere_names:
    if t_name != "Y1p":
        short_distances[t_name].sort(reverse=True, key=key_for_sort)
        minimum_distance = short_distances[t_name][2][1]
        long_distances[t_name].sort(key=key_for_sort)
        maximum_distance = long_distances[t_name][2][1]
        for scankey in scans:
            scan = scans[scankey]
            t_mere = scan.get(t_name)
            if t_mere != None:
                if t_mere["dist2"] > maximum_distance:
                    t_mere["z_sign"] = "positive"
                elif t_mere["dist1"] < minimum_distance:
                    t_mere["z_sign"] = "negative"
                else:
                    t_mere["z_sign"] = "ambiguous"

print(scans)

              

            
