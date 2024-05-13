import sys
import csv
import math
from copy import deepcopy

scans = {} # empty dictionary.  Will contain one entry per scan.
telomere_names=[]
skip_these = ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]
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
                    if not t_name in skip_these:
                        scan_values[t_name] = t_values # will contain all descriptors of this telomere
                        t_values["x"] = float(str)
                elif row_pos>1 and (row_pos-2) % 5 == 0 and not t_name in skip_these:
                    t_values["y"] = float(str)
                elif row_pos>2 and (row_pos-3) % 5 == 0 and not t_name in skip_these:
                    t_values["real_z"] = float(str)
                row_pos += 1
        row_count += 1
#  print(scans)

for t_name in ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]:
    if t_name in telomere_names:
        telomere_names.remove(t_name)

# clean the data to remove any outliers

for scankey in scans:
    scan = scans[scankey]
    bad_points = []
    for t_name in scan:
        t_mere = scan[t_name]
        if t_mere["x"]**2 + t_mere["y"]**2 > 25:
            len = t_mere["x"]**2 + t_mere["y"]**2
            x = t_mere["x"]
            y = t_mere["y"]
            print(f"found bad point in scan {scankey} {t_name} {x} {y} {len} ")
            bad_points.append(t_name)
    for t_name in bad_points:
        scan.pop(t_name)
            

# we choose a reference point.  This will be the one nearest the center, as that
# best disambiguates the rest

def key_for_sort(e):
    return e[1]

def key_for_sort2(e):
    return e[2]



# Assign z values -- all positive for the moment

for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        t_mere = scan[t_name]
        t_mere["z"] = math.sqrt(25 - t_mere["x"]**2 - t_mere["y"]**2)
        if abs(t_mere["z"]-abs(t_mere["real_z"])) >0.1:
            print(t_mere)
            raise Exception("big error here!")

# We will assume that Y1p has a z value that is positive.  This might be wrong, or at least,
# inconsistent with the others, in which case we reflect the scan at the end.

# here are some random choices we might adjust
refpoint1 = "Y1p"
pts_to_agree = 1
experimental_error = 0.5
reflection = 1

for scankey in scans:
    scan = scans[scankey]
    scan[refpoint1]["z_sign"] = "positive"  # might be false -- in case we reflect at the end
    scan[refpoint1]["arbitrary"] = True
    if scan[refpoint1]["real_z"] < 0:
        reflection = -1
    for t_name in scan:
        if t_name != refpoint1:
            scan[t_name]["z_sign"] = "ambiguous" # all others ambiguous at the start

count = 0
ambiguous = 0
def sign_of(n):
    if n>=0:
        return True
    return False
def eliminate_ambiguous(refpoint):
    global count
    global ambiguous
    # Ok, choose the other values for z_sign that best match.  See logic below
    distances = {}
    for scankey in scans:
        scan = scans[scankey]
        if scan.get(refpoint) == None:
            continue
        refp = deepcopy(scan[refpoint])
        if refp["z_sign"] == "ambiguous": # unlikely to resolve anything
            # unless all are ambiguous except for the initial one
            all_ambiguous = True
            arbitrary = refpoint1
            for t_name in scan:
                if scan[t_name].get("arbitrary") != None:
                    arbitrary = t_name
                elif scan[t_name]["z_sign"] != "ambiguous":
                    all_ambiguous = False
                    break
            if all_ambiguous:
                scan[refpoint]["z_sign"] = "positive"
                scan[refpoint]["arbitrary"] = True
                scan[arbitrary]["z_sign"] = "ambiguous"
                scan[arbitrary].pop("arbitrary")
            else:
                continue    
        if refp["z_sign"] == "negative":  # I think this is not right
            continue # for now 

        for t_name in scan:
           # z_sign = scan[t_name]["z_sign"]   
            if t_name != refpoint and scan[t_name]["z_sign"] == "ambiguous":
                if distances.get(t_name) == None:
                    distances[t_name] = []
                t_mere = scan[t_name]
                distance1 = (refp["x"] - t_mere["x"])**2 + (refp["y"] - t_mere["y"])**2 + (refp["z"] - t_mere["z"])**2
                distance2 = (refp["x"] - t_mere["x"])**2 + (refp["y"] - t_mere["y"])**2 + (refp["z"] + t_mere["z"])**2
                # if refp["z_sign"] == "negative": # the first values will be longer in this case
                #     distance_save = distance2
                #     distance2 = distance1
                #     distance1 = distance_save
                distance3 = distance2 - distance1
                distances[t_name].append((scankey, distance3, distance1, distance2)) 
                if scankey == "test20" and t_name== "Xp":
                    print(distance1)
                    print(distance2)
                    print(t_mere["z"])
                    print(refp["z"])

    for t_name in distances:
        if (t_name=="Xp"):
            for distance in distances[t_name]:
                if distance[0] == "test20":
                    print(distance)
        distances[t_name].sort(key = key_for_sort) 
        if (t_name=="Xp"):
            for distance in distances[t_name]:
                if distance[0] == "test20":
                    print(distance)
        # The points at the top of this sorted list will be the ones closest to the edge of
        # the circle, and therefore the ones for which the computed distance will vary
        # the least depending on the sign of the z_value.  It is unlikely that distances
        # smaller than experimental error less than the lesser of these two distances
        # are correct.  It is unlikely that distances greater than experimental error
        # plus the greater of these two distances are correct.  So, we can eliminate
        # some distances, but many remain ambiguous
        found_error = False
        compare_distance = distances[t_name][pts_to_agree-1]
        for distance in distances[t_name][pts_to_agree:]:
            if distance[3] > compare_distance[3] and distance[2] < compare_distance[2]:
                continue # not diagnostic, we can't tell
            if distance[3] > compare_distance[3]  + experimental_error:
                scans[distance[0]][t_name]["z_sign"] = "positive"
                if scans[distance[0]][t_name]["real_z"]*reflection < 0:
                    found_error = True
            elif distance[2] < distances[t_name][pts_to_agree-1][2] - experimental_error:
                scans[distance[0]][t_name]["z_sign"] = "negative"
                if scans[distance[0]][t_name]["real_z"]*reflection > 0:
                    found_error = True
            if found_error and not scans[compare_distance[0]][t_name]["z_sign"]== "ambiguous":
                print("starting error report")
                print(f"t_name: {t_name}")
                print(f"distance in question: {distance}")
                print(f"distance to compare: {compare_distance}")
                print(reflection)
                find_point = scans[distance[0]][t_name]
                find_ref = scans[distance[0]][refpoint]
                compare_point = scans[compare_distance[0]][t_name]
                compare_ref = scans[compare_distance[0]][refpoint]
                
                print(f"point to find: {scans[distance[0]][t_name]}")
                print(f"ref point: {scans[distance[0]][refpoint]}")
                print(f"compare point: {scans[compare_distance[0]][t_name]}")
                print(f"compare ref: {scans[compare_distance[0]][refpoint]}")
                find_z = -abs(find_point["z"])
                if sign_of(find_point["real_z"]):
                    find_z = abs(find_point["z"])
                find_rz = -abs(find_ref["z"])        
                if sign_of(find_ref["real_z"]):
                    find_rz = abs(find_ref["z"])
                compare_z = -abs(compare_point["z"])
                if sign_of(compare_point["real_z"]):
                    compare_z = abs(compare_point["z"])
                compare_rz = -abs(compare_ref["z"])
                if sign_of(compare_ref["real_z"]):
                    compare_rz = abs(compare_ref["z"])
                fprd = (find_point["x"]-find_ref["x"])**2 + (find_point["y"] - find_ref["y"])**2 \
                    + (find_z - find_rz)**2
                cprd =  (compare_point["x"]-compare_ref["x"])**2 + (compare_point["y"] - compare_ref["y"])**2 \
                    + (compare_z - compare_rz)**2               
                print(f"actual distances for find and compare points: {fprd} {cprd}")
                raise Exception("Something went wrong!")
                
    count = 0
    ambiguous = 0
    for scankey in scans:
        scan = scans[scankey]
        for t_name in scan:
            count += 1
            if scan[t_name]["z_sign"] == "ambiguous":
                z = scan[t_name]["z"]
                z2 = scan[refpoint]["z"]
                print(f"ambiguous {scankey} {refpoint} {t_name} {z} {z2}")
                ambiguous += 1
    print(f"Of {count} z values, {ambiguous} are still ambiguous.")

eliminate_ambiguous(refpoint1)

pts_to_try = []

for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        if t_name != refpoint1 and not t_name in pts_to_try:
            pts_to_try.append(t_name)

while ambiguous > 0 and len(pts_to_try) > 0:
    eliminate_ambiguous(pts_to_try[0])
    pts_to_try.pop(0)

# That's about all we can tell
        
# Now make tuples for each point

for scankey in scans:
    scan = scans[scankey]
    for t_name in scan:
        t_mere = scan[t_name]
        if t_mere["z_sign"] == "positive":
            t_mere["loc"] = (t_mere["x"],t_mere["y"],t_mere["z"])
        elif t_mere["z_sign"] == "negative":
            t_mere["loc"] = (t_mere["x"],t_mere["y"],-t_mere["z"])
        else: 
            t_mere["loc"] = (t_mere["x"],t_mere["y"],t_mere["z"])
            t_mere["loc_alt"] = (t_mere["x"],t_mere["y"],-t_mere["z"])
                     
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
    if y1p[0] == 0 and y1p[1] > 0:
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
    if scan["Y1p"]["loc"][1]<0:
        print(f"problem here in loc, scan is {scankey}")
    matrix3 = first_matrix(scan["Xp"]["loc"])
    for t_name in scan:
        t_mere = scan[t_name]
        t_mere["loc"] = matrix_rotate(matrix3,t_mere["loc"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == matrix_rotate(matrix3, t_mere["loc_alt"])


# ok, now we need to see if reflection is necessary.  We need to find the telomere
# for which the absolute value of the z is largest

# z_abs = {}
# for scankey in scans:
#     scan = scans[scankey]
#     for t_name in scan:
#         t_mere = scan[t_name]
#         if t_mere["z_sign"] != "ambiguous":  
#             if z_abs.get(t_name) == None:
#                 z_abs[t_name] = [0,0]
#             z_abs[t_name][0] += 1
#             z_abs[t_name][1] += abs(t_mere["loc"][2])

# z_avg=[]
# for t_name in z_abs:
#     z_ab = z_abs[t_name]
#     z_avg.append((t_name,z_ab[1]/float(z_ab[0])))

# z_avg.sort(reverse = True, key = key_for_sort)
# z_max_t_name = z_avg[0][0]

# # having found the largest average z, we now make sure that that is positive in each
# # scan. If we have to change the sign, we change the sign for all the z values in the scan.
# # This performs the reflection.



# bad_scans = []
# for scankey in scans:
#     scan = scans[scankey]
#     t_mere = scan.get(z_max_t_name)
#     if t_mere == None:
#         # hm.  This is a problem.  Of course this shouldn't happen
#         # but if it does we can fix it by reference to some other telomere
#         print("found a bad scan")
#         bad_scans.append(scankey)
#     elif t_mere["z_sign"] == "ambiguous" and sign_of(t_mere["loc"][2]) != sign_of(t_mere["loc_alt"][2]):
#         print("found another bad scan")
#         bad_scans.append(scankey)
#     elif not sign_of(t_mere["loc"][2]): # if the telomere z is negative, we have to reflect it
#         print("reflecting")
#         for t_name in scan:
#             t_mere = scan[t_name]
#             new_loc = (t_mere["loc"][0],t_mere["loc"][1],-t_mere["loc"][2])
#             t_mere["loc"] = new_loc
#             if t_mere["z_sign"] == "ambiguous": # could happen, if a rotation put both possible locations here
#                 new_loc = (t_mere["loc_alt"][0],t_mere["loc_alt"][1],-t_mere["loc_alt"][2])
#                 t_mere["loc_alt"] = new_loc

# # because of data errors, bad_scans will not be empty.  So we need to pick another likely value, and
# # figure out the sign for it.  We can't just set it to positive -- instead we set it to the 
# # sign that it has in the majority of cases already resolved.
# if len(bad_scans) > 0:
#     sign_count = 0
#     for scankey in scans:
#         if not scankey in bad_scans:
#             scan = scans[scankey]
#             t_mere = scan[z_avg[2][0]]
#             if t_mere["z_sign"] == "positive":
#                 sign_count += 1
#             elif t_mere["z_sign"] == "negative":
#                 sign_count -=1

#     # hopefully sign_count won't be 0!
#     if sign_count == 0:
#         print("gee, how could that happen?")

#     if sign_count > 0:
#         for scankey in bad_scans:
#             scan = scans[scankey]
#             if scan[z_avg[2][0]]["z_sign"] == "ambiguous":
#                 print("man, we ain't got no luck.")
#             elif scan[z_avg[2][0]]["z_sign"] == "negative": #gotta flip it
#                 for t_name in scan:
#                     t_mere = scan[t_name]
#                     new_loc = (t_mere["loc"][0], t_mere["loc"][1], -t_mere["loc"][2])
#                     t_mere["loc"] = new_loc
#     else:
#         for scankey in bad_scans:
#             scan = scans[scankey]
#             if scan[z_avg[2][0]]["z_sign"] == "ambiguous":
#                 print("man, we ain't got no luck.")
#             elif scan[z_avg[2][0]]["z_sign"] == "positive": #gotta flip it
#                 for t_name in scan:
#                     t_mere = scan[t_name]
#                     new_loc = (t_mere["loc"][0], t_mere["loc"][1], -t_mere["loc"][2])
#                     t_mere["loc"] = new_loc
        
# to figure out reflections, we orient each scan so as to put the majority of points in
# the same hemispheres

t_names = ["Y1p", "Xp", "Xq", "Y1q", "Y2p", "Y2q"]
t_counts = {}
for t_name in t_names:
    for scankey in scans:
        scan = scans[scankey]
        if scan.get(t_name) != None:
            if scan[t_name]["z_sign"] != "ambiguous" \
                     and abs(scan[t_name]["z"])>1.0:
                if t_counts.get(t_name) == None:
                    t_counts[t_name] = 0
                if scan[t_name]["loc"][2] > 0:
                    t_counts[t_name] += 1
                else:
                    t_counts[t_name] -= 1

print(t_counts)

to_reflect = []
ambiguous = []

def reflect(scan):


for scankey in scans:
    scan = scans[scankey]
    if scan.get('Xp') != None and scan["Xp"]["loc"][2] > 1:


raise Exception("stop here")

for scankey in scans:
    scan = scans[scankey]
    #if scan.get("Y1p") != None:

    for t_name in scan:
            t_mere = scan[t_name]
            new_loc = (t_mere["loc"][0],t_mere["loc"][1],-t_mere["loc"][2])
            t_mere["loc"] = new_loc
            if t_mere["z_sign"] == "ambiguous": # could happen, if a rotation put both possible locations here
                new_loc = (t_mere["loc_alt"][0],t_mere["loc_alt"][1],-t_mere["loc_alt"][2])
                t_mere["loc_alt"] = new_loc
raise Exception("stop here")



        
        
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
    d = std_deviations[t_name]["stddev"]
    print(f"{t_name} {d}")

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

    


              

            
