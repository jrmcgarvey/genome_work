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
    ref_point["z_sign"] = "positive"
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

# print(scans)

# try once again to disambiguate any ambiguous Xp

bad_xps = []
for scankey in scans:
    if scans[scankey]["Xp"]["z_sign"] == "ambiguous":
        bad_xps.append(scankey) # this list will have all the scans where Xp is ambiguous

if len(bad_xps) > 0: # still some to disambiguate, this time w/r reference to Y1q
    short_lengths = []
    long_lengths = []
    for scankey in scans:
        scan = scans[scankey]
        xp = scan["Xp"]
        y1q = scan.get("Y1q")
        if y1q != None:
            y1q = scan["Y1q"]
            short_length = (xp["z1"] - y1q["z1"])**2 \
                     + (xp["x"] - y1q["x"])**2 \
                     + (xp["y"] - y1q["y"])**2 
            short_lengths.append(short_length)
            long_length = (xp["z1"] - y1q["z2"])**2 \
                     + (xp["x"] - y1q["x"])**2 \
                     + (xp["y"] - y1q["y"])**2 
            long_lengths.append(long_length)
    short_lengths.sort(reverse=True)
    long_lengths.sort()
    max_len = long_lengths[2]
    min_len = short_lengths[2]
    for scankey in bad_xps:
        scan = scans[scankey]
        xp = scan["Xp"]
        y1q = scan["Y1q"]
        short_value = (xp["z1"] - y1q["z1"])**2 \
                     + (xp["x"] - y1q["x"])**2 \
                     + (xp["y"] - y1q["y"])**2 
        long_value = (xp["z2"] - y1q["z1"])**2 \
                     + (xp["x"] - y1q["x"])**2 \
                     + (xp["y"] - y1q["y"])**2 
        if long_value > max_len:
            scans[scankey]["Xp"]["z_sign"] = scans[scankey]["Y1q"]["z_sign"]
        elif short_value < min_len:
            if scans[scankey]["Y1q"]["z_sign"] == "positive":
                scans[scankey]["Xp"]["z_sign"] = "negative"
            else:
                scans[scankey]["Xp"]["z_sign"] = "positive"
        else:
            print("Discarding " + scankey + " as Xp z value is ambiguous")
            scans.pop(scankey)
        
# Now make tuples for each point

for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        t_mere = scan[t_name]
        if t_mere["z_sign"] == "positive":
            t_mere["loc"] = (t_mere["x"],t_mere["y"],t_mere["z1"])
        elif t_mere["z_sign"] == "negative":
            t_mere["loc"] = (t_mere["x"],t_mere["y"],t_mere["z2"])
        else: 
            t_mere["loc"] = (t_mere["x"],t_mere["y"],t_mere["z1"])
            t_mere["loc_alt"] = (t_mere["x"],t_mere["y"],t_mere["z2"])
                     
# Now a function to multiply a 3x3 matrix times a tuple

def matrix_rotate(m,t):
    x = m[0][0] * t[0] + m[0][1] * t[1] + m[0][2] * t[2]
    y = m[1][0] * t[0] + m[1][1] * t[1] + m[1][2] * t[2]
    z = m[2][0] * t[0] + m[2][1] * t[1] + m[2][2] * t[2]
    return (x,y,z)

# Now a function to calculate the first matrix from the Y1p position.  This
# rotation is about the Y axis, and moves Y1p to a point in the x, y plane where x is positive.
# The rotation is the identity if Y is already in that plane with a non-negative x value

def first_matrix(y1p):
    if y1p[2] == 0:
        if y1p[0] >= 0:
            return [[1,0,0],[0,1,0],[0,0,1]]
        else:
            return [[-1,0,0],[0,1,0],[0,0,-1]]
    else:
        hyp = math.sqrt(y1p[2]**2 + y1p[0]**2)
        sin_theta = y1p[2] / hyp
        cos_theta = y1p[0] / hyp
        return [[cos_theta, 0, sin_theta], [0,1,0],[-sin_theta, 0, cos_theta]]

# Next we want to rotate about the z axis, to bring Y1p to the north pole

def second_matrix(y1p):
    if (y1p[0]) == 0:
        return [[1,0,0],[0,1,0],[0,0,1]]
    else:
        sin_theta = y1p[0] / 5 # 5 being the radius
        cos_theta = y1p[1] / 5
        return [[cos_theta, -sin_theta, 0], [sin_theta, cos_theta, 0], [0,0,1]]

# The third matrix is just like the first, a rotation about the y axis, but this time
# computed using Xp.  So we don't need another function.
#diagnostics, to be removed
# val = (1,1,1)

# def standard_length(tup):
#     ra = []
#     for num in tup:
#         ra.append(num * 5 / math.sqrt(tup[0]**2 + tup[1]**2 + tup[2]**2))
#     print(ra)
#     return ra

# val2 = standard_length(val)

# m1 = first_matrix(val2)
# val3 = matrix_rotate(m1,val2)
# print(val3)

# m2 = second_matrix(val3)

# val4 = matrix_rotate(m2, val3)
# print(val4)
# Ok, now we compute and apply each matrix for each scan:

for scankey in scans:
    scan = scans[scankey]
    matrix1 = first_matrix(scan["Y1p"]["loc"])
    for t_name in scan:
        t_mere = scan[t_name]
        t_mere["loc"] = matrix_rotate(matrix1,t_mere["loc"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == matrix_rotate(matrix1,t_mere["loc_alt"])
    matrix2 = second_matrix(scan["Y1p"]["loc"])
    for t_name in scan:
        t_mere = scan[t_name]
        t_mere["loc"] = matrix_rotate(matrix2,t_mere["loc"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == matrix_rotate(matrix2,t_mere["loc_alt"])
    matrix3 = first_matrix(scan["Xp"]["loc"])
    for t_name in scan:
        t_mere = scan[t_name]
        t_mere["loc"] = matrix_rotate(matrix3,t_mere["loc"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == matrix_rotate(matrix3, t_mere["loc_alt"])
        print("values for scan " + scankey)
        print(scan["Y1p"]["loc"])
        print(scan["Xp"]["loc"])

# ok, now we need to see if reflection is necessary.  We need to find the telomere
# for which the absolute value of the z is largest

z_abs = {}
for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        t_mere = scan[t_name]
        if t_mere["z_sign"] != "ambiguous":  
            if z_abs.get(t_name) == None:
                z_abs[t_name] = [0,0]
            z_abs[t_name][0] += 1
            z_abs[t_name][1] += abs(t_mere["loc"][2])

z_avg=[]
for t_name in z_abs:
    z_ab = z_abs[t_name]
    z_avg.append((t_name,z_ab[1]/float(z_ab[0])))

z_avg.sort(reverse = True, key = key_for_sort)
z_max_t_name = z_avg[0][0]

# having found the largest average z, we now make sure that that is positive in each
# scan. If we have to change the sign, we change the sign for all the z values in the scan.
# This performs the reflection.

def sign_of(n):
    if n>=0:
        return True
    return False

bad_scans = []
for scankey in scans:
    scan = scans[scankey]
    t_mere = scan.get(z_max_t_name)
    if t_mere == None:
        # hm.  This is a problem.  Of course this shouldn't happen
        # but if it does we can fix it by reference to some other telomere
        print("found a bad scan")
        bad_scans.append(scankey)
    elif t_mere["z_sign"] == "ambiguous" and sign_of(t_mere["loc"][2]) != sign_of(t_mere["alt_loc"][2]):
        print("found another bad scan")
        bad_scans.append(scankey)
    elif not sign_of(t_mere["loc"][2]): # if the telomere z is negative, we have to reflect it
        for t_name in scan:
            t_mere = scan[t_name]
            new_loc = (t_mere["loc"][0],t_mere["loc"][1],-t_mere["loc"][2])
            t_mere["loc"] = new_loc
            if t_mere["z_sign"] == "ambiguous": # could happen, if a rotation put both possible locations here
                new_loc = (t_mere["loc_alt"][0],t_mere["loc_alt"][1],-t_mere["loc_alt"][2])
                t_mere["loc_alt"] = new_loc

# because of data errors, bad_scans will not be empty.  So we need to pick another likely value, and
# figure out the sign for it.  We can't just set it to positive -- instead we set it to the 
# sign that it has in the majority of cases already resolved.

sign_count = 0
for scankey in scans:
    if not scankey in bad_scans:
        scan = scans[scankey]
        t_mere = scan[z_avg[2][0]]
        if t_mere["z_sign"] == "positive":
            sign_count += 1
        elif t_mere["z_sign"] == "negative":
            sign_count -=1

# hopefully sign_count won't be 0!
if sign_count == 0:
    print("gee, how could that happen?")

if sign_count > 0:
    for scankey in bad_scans:
        scan = scans[scankey]
        if scan[z_avg[2][0]]["z_sign"] == "ambiguous":
            print("man, we ain't got no luck.")
        elif scan[z_avg[2][0]]["z_sign"] == "negative": #gotta flip it
            for t_name in scan:
                t_mere = scan[t_name]
                new_loc = (t_mere["loc"][0], t_mere["loc"][1], -t_mere["loc"][2])
                t_mere["loc"] = new_loc
else:
    for scankey in bad_scans:
        scan = scans[scankey]
        if scan[z_avg[2][0]]["z_sign"] == "ambiguous":
            print("man, we ain't got no luck.")
        elif scan[z_avg[2][0]]["z_sign"] == "positive": #gotta flip it
            for t_name in scan:
                t_mere = scan[t_name]
                new_loc = (t_mere["loc"][0], t_mere["loc"][1], -t_mere["loc"][2])
                t_mere["loc"] = new_loc
        
        
# ok, that's it for the reflections -- provided we don't find bad scans.
# Now for the final step -- best fit. We're going to have to find the average location twice,
# so I'll create a function

def find_center(scans):
    t_mere_locs = {}
    for scankey in scans:
        scan = scans[scankey]
        for t_name in scan:
            t_mere = scan[t_name]
            if t_mere["z_sign"] != "ambiguous":
                if t_mere_locs.get(t_name) == None:
                    t_mere_locs[t_name] = [t_mere["loc"]]
                else: 
                    t_mere_locs[t_name].append(t_mere["loc"])
    t_mere_centers = {}
    for t_name in t_mere_locs:
        x_tot = 0
        y_tot = 0
        z_tot = 0
        loc_array = t_mere_locs[t_name]
        for loc in loc_array:
            x_tot += loc[0]
            y_tot += loc[1]
            z_tot += loc[2]
        x_mean = x_tot / float(len(loc_array))
        y_mean = y_tot / float(len(loc_array))
        z_mean = z_tot / float(len(loc_array))
        multiplier = 5 / math.sqrt(x_mean**2 + y_mean**2 + z_mean**2)
        t_mere_centers[t_name] = (x_mean * multiplier, y_mean * multiplier, z_mean * multiplier)
        # this last is to put it on the surface of the sphere.
    return t_mere_centers

t_mere_centers = find_center(scans) # first time
print(t_mere_centers)

# Now, best fit for any points left ambiguous
for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        t_mere = scan[t_name]
        if t_mere["z_sign"] == "ambiguous": # find the one that's closest
            center = t_mere_centers[t_name]
            dist1 = (center[0] - t_mere["loc"][0]) ** 2 + (center[1] - t_mere["loc"][1]) ** 2 \
                  + (center[2] - t_mere["loc"][2]) ** 2
            dist2 = (center[0] - t_mere["loc_alt"][0]) ** 2 + (center[1] - t_mere["loc_alt"][1]) ** 2 \
                  + (center[2] - t_mere["loc_alt"][2]) ** 2
            if dist1 > dist2:
                t_mere["loc"] = t_mere["loc_alt"]
            t_mere["z_sign"] = "fixed"

# recompute the centers to match

t_mere_centers = find_center(scans)  

# ok, now to find the standard deviations

std_deviations = {}

for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        if std_deviations.get(t_name) == None:
            std_deviations[t_name] = {"sum": float(0), "count": 0}
        center = t_mere_centers[t_name]
        loc = scan[t_name]["loc"]
        std_deviations[t_name]["sum"] += (center[0] - loc[0])**2 + (center[1]-loc[1])**2 + (center[2]-loc[2])**2 
        std_deviations[t_name]["count"] += 1

for t_name in std_deviations:
    std_deviations[t_name]["stddev"] = math.sqrt(std_deviations[t_name]["sum"]/float(std_deviations[t_name]["count"]))

print(std_deviations)

f = open(sys.argv[1]+".result", "w")

for t_name in t_mere_centers:
    center = t_mere_centers[t_name]
    stddev = std_deviations[t_name]
    x = center[0]
    y = center[1]
    z = center[2]
    w = std_deviations[t_name]["stddev"]
    f.write(f"{t_name} x: {x:.2f} y: {y:.2f} z: {z:.2f} stddev: {w:.2f}\n")
f.close()
    
    
    


              

            
