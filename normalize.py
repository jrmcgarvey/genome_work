import sys
import csv
import math
from copy import deepcopy
import utils

scans = {} # empty dictionary.  Will contain one entry per scan.
telomere_names=[]
skip_these = ["1Ap", "1Aq", "2Ap", "2Aq", "1Bp", "1Bq", "2Bp", "2Bq"]
includez = False
if sys.argv[2] == "includez": # for testing
    includez = True
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
            scan = {}
            scans[row[0]] = scan # will contain a dictionary of telomeres
            row_pos = 0
            for str in row:
                t_name = ""
                if row_pos>0 and (row_pos-1) % 5 == 0:
                    t_name = telomere_names[int((row_pos-1)/5)]
                    if not t_name in skip_these:
                        t_mere={}
                        scan[t_name]=t_mere
                        p=[0,0,0]
                        t_mere["p"]=p
                        p[0] = float(str)
                elif row_pos>1 and (row_pos-2) % 5 == 0 and not t_name in skip_these:
                    p[1] = float(str)
                elif includez and row_pos>2 and (row_pos-3) % 5 == 0 and not t_name in skip_these:
                     t_mere["real_z"]= float(str)
                row_pos += 1
        row_count += 1


for t_name in skip_these:
    if t_name in telomere_names:
        telomere_names.remove(t_name)

# clean the data to remove any outliers

for scankey, scan in scans.items():
    bad_points = []
    for t_name in scan:
        t_mere = scan[t_name]
        p=t_mere["p"]
        len2 =  p[0]**2 + p[1]**2
        if len2>25.5: # give the measurement a little slack
            print(f"found bad point in scan {scankey} {t_name} {p[0]} {p[1]} {len} ")
            bad_points.append(t_name)
        elif len2>25:
            p[0] = p[0] / math.sqrt(25/len2)
            p[1] = p[1] / math.sqrt(25/len2)
            p[2] = 0
        else:
            p[2] = math.sqrt(25 - p[0]**2 - p[1]**2)
        if abs(p[2] - abs(t_mere["real_z"])) > 0.1:
            raise Exception("bad error in test {t_mere}")
    for t_name in bad_points:
        scan.pop(t_name)
            

def key_for_sort(e):
    return e[1]

def key_for_sort2(e):
    return e[2]

def reflect(scan):
    for t_mere in scan.values():
        if t_mere["z_sign"] == "positive":
            t_mere["z_sign"] = "negative"
        elif t_mere["z_sign"] == "negative":
            t_mere["z_sign"] = "positive"

pts_to_agree = 1
experimental_error = 0.1
count = 0 
ambiguous = 0
def sign_of(n):
    if n>=0:
        return 1
    return -1
def eliminate_ambiguous(refpoint):
    global count
    global ambiguous
    incorrect = 0
    arbitrary = []
    # Ok, choose the other values for z_sign that best match.  See logic below
    for scan in scans.values():
        if scan.get(refpoint) == None:
            continue
        t_mere = scan[refpoint]
        if t_mere["z_sign"] == "ambiguous":
            t_mere["z_sign"] = "positive" #tentative, for a start
            t_mere["arbitrary"] = True
            arbitrary.append(t_mere)
        elif t_mere["z_sign"] == "negative": # resolved this way from another point
            reflect(scan) # we need to calculate this way
    distances = {}
    for scankey, scan in scans.items():
        if scan.get(refpoint) == None:
            continue
        refp = scan[refpoint]
        if refp["p"][2]<experimental_error:
            continue
        for t_name,t_mere in scan.items():
            if t_mere["p"][2] < experimental_error:
                continue
           # z_sign = scan[t_name]["z_sign"]   
            if t_name != refpoint and t_mere["z_sign"] == "ambiguous":
                if distances.get(t_name) == None:
                    distances[t_name] = []
                p_ref=deepcopy(refp["p"])
                p_c = t_mere["p"]
                if refp["z_sign"]: 
                    p_ref[2] *= -1
                distance1 = (p_ref[0] - p_c[0])**2 + (p_ref[1] - p_c[1])**2  + (p_ref[2] - p_c[2])**2
                distance2 = (p_ref[0] - p_c[0])**2 + (p_ref[1] - p_c[1])**2  + (p_ref[2] + p_c[2])**2
                # if refp["z_sign"] == "negative": # the first values will be longer in this case
                #     distance_save = distance2
                #     distance2 = distance1
                #     distance1 = distance_save
                distance3 = distance2 - distance1
                distances[t_name].append((scankey, distance3, distance1, distance2))

    for t_name in distances:
        distances[t_name].sort(key = key_for_sort) 
        for distance in distances[t_name]:

        # The points at the top of this sorted list will be the ones closest to the edge of
        # the circle, and therefore the ones for which the computed distance will vary
        # the least depending on the sign of the z_value.  It is unlikely that distances
        # smaller than experimental error less than the lesser of these two distances
        # are correct.  It is unlikely that distances greater than experimental error
        # plus the greater of these two distances are correct.  So, we can eliminate
        # some distances, but many remain ambiguous
            compare_distance = distances[t_name][pts_to_agree-1]
            for distance in distances[t_name][pts_to_agree:]:
                scankey = distance[0]
                scan = scans[scankey]
                t_mere = scan[t_name]
                t_mere_ref = scan[refpoint]
                compare_ref = scans[compare_distance[0]][refpoint]
                if abs(distance[3] - distance[2]) < experimental_error:
                    continue # can't resolve
                if distance[3] > compare_distance[3]- experimental_error/2 and distance[2] < compare_distance[2]+experimental_error/2:
                    continue # not diagnostic, we can't tell
                if distance[3] < compare_distance[3] and distance[2] > compare_distance[2]:
                    continue # also not diagnostic
                if distance[3] > compare_distance[3]  + experimental_error:
                    t_mere["z_sign"] = "positive"
                elif distance[2] < distances[t_name][pts_to_agree-1][2] - experimental_error:
                    t_mere["z_sign"] = "negative"
                if includez and t_mere["z_sign"] != "ambiguous":
                    sign1 = 1
                    if t_mere["z_sign"] == "negative":
                        sign1 = -1
                    sign2 = sign_of(t_mere["real_z"])
                    sign3 = sign_of(t_mere_ref["real_z"])
                    sign4 = 1
                    if t_mere_ref["z_sign"]=="negative":
                        sign4 = -1
                    if (sign1 * sign2 * sign3 * sign4) <0:
                        # this suggests an error.  If z_sign is negative, we mean that
                        # the point is not in the same z hemisphere as the reference point for
                        # that scan.  That should be true only if the real_z signs of the 
                        # point and the reference point differ.  
                        print("starting error report")
                        print(f"t_name: {t_name}")
                        print(f"distance in question: {distance}")
                        print(f"distance to compare: {compare_distance}")
                        t_mere_ref = scan[refpoint]
                        scankey_c = compare_distance[0]
                        scan_r = scans[scankey_c]
                        t_mere_c = scan_r[t_name]
                        t_mere_c_ref = scan_r[refpoint]
                        print(f"t_mere: {t_mere}")
                        print(f"t_mere_ref: {t_mere_ref}")
                        print(f"t_mere_c: {t_mere_c}")
                        print(f"t_mere_c_ref: {t_mere_c_ref}")
                        p = t_mere["p"]
                        z1 = t_mere["real_z"]
                        p_r = t_mere_ref["p"]
                        z_r = t_mere_ref["real_z"]
                        dist1 = (p[0] - p_r[0])**2 + (p[1]-p_r[1])**2 + (z1 -z_r)**2
                        p = t_mere_c["p"]
                        z1 = t_mere_c["real_z"]
                        p_r = t_mere_c_ref["p"]
                        z_r = t_mere_c_ref["real_z"]
                        dist2 = (p[0] - p_r[0])**2 + (p[1]-p_r[1])**2 + (z1 -z_r)**2
                        print(f"reported distances (squared) are: {dist1} {dist2}")
                        raise Exception("Something went wrong!")
    for t_mere in arbitrary:
        t_mere["z_sign"] = "ambiguous"
        t_mere.pop("arbitrary")                
    count = 0
    ambiguous = 0
    for scankey, scan in scans.items():
        for t_name,t_mere in scan.items():
            count += 1
            if  t_mere["z_sign"] == "ambiguous":
                z = t_mere["p"][2]
                z2 = scan[refpoint]["p"][2]
                ambiguous += 1
            elif includez:
                sign2 = 1
                if t_mere["z_sign"] == "negative":
                    sign2 = -1
                sign3 = sign_of(t_mere["real_z"])
                if sign2 * sign3 < 0:
                    incorrect += 1
    print(f"Of {count} z values, {ambiguous} are still ambiguous.")


for scan in scans.values():
    for t_mere in scan.values():
        t_mere["z_sign"] = "ambiguous" # all ambiguous at the start

for t_name in telomere_names:
    eliminate_ambiguous(t_name)


# pts_to_try = []

# for scan in scans.values():
#     for t_name in scan:
#         if t_name != refpoint1 and not t_name in pts_to_try:
#             pts_to_try.append(t_name)

# while ambiguous > 0 and len(pts_to_try) > 0:
#     eliminate_ambiguous(pts_to_try[0])
#     pts_to_try.pop(0)

# That's about all we can tell
        
# Now make tuples for each point

for scan in scans.values():
    for t_mere in scan.values():
        if t_mere["z_sign"] == "negative":
            t_mere["p"][2] *= -1
        elif t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] = deepcopy(t_mere["p"])
            t_mere["loc_alt"][2] *= -1

# for scan in scans.values():
#     for t_mere in scan.values():
#         t_mere["p"][2] *= sign_of(t_mere["real_z"])
#         p_s = "p"
#         real_z = "real_z"
#         #print(f"z values {t_mere[p_s][2]} {t_mere[real_z]}")
#         #t_mere["p"][2] = t_mere["real_z"]
#         t_mere["z_sign"] = "fixed"
#raise Exception("quit here")
bad_scans=[]
north_pole = telomere_names[0]
xy_plane = telomere_names[1]
for scankey, scan in scans.items():
    if scan.get(north_pole) == None:
        print("no north pole")
        bad_scans.append(scankey)
        continue
    if scan.get(xy_plane) == None:
        print("no xy_plane")
        bad_scans.append(scankey)
        continue
    if scan[xy_plane]["z_sign"] == "ambiguous":
        print("xy_plane ambiguous")
        bad_scans.append(scankey)
        continue
    if scan[north_pole]["z_sign"] == "ambiguous": # doesn't happen, but we might change the logic
        print("north pole ambiguous")
        bad_scans.append(scankey)
        continue # we might be able to put these back in, by rotating w/r to different points
    # but for now we will just toss these.
    matrix1 = utils.rotate_to_x_y(scan[north_pole]["p"])
    p_s = "p"
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix1,t_mere["p"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix1,t_mere["loc_alt"])
    print(f"scan np: {scan[north_pole][p_s]}")
    matrix2 = utils.rotate_to_y_z(scan[north_pole]["p"])
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix2,t_mere["p"])

        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix2,t_mere["loc_alt"])
    print(f"scan np: {scan[north_pole][p_s]}")
    print(f"scan xy: {scan[xy_plane][p_s]}")
    matrix3 = utils.rotate_to_x_y(scan[xy_plane]["p"])
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix3,t_mere["p"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix3, t_mere["loc_alt"])
    print(matrix3)
    print(f"scan np: {scan[north_pole][p_s]}")
    print(f"scan xy: {scan[xy_plane][p_s]}")
#    raise Exception("stop here")
print(f"threw away {len(bad_scans)} scans")
bad_scans_save = {}
for scankey in bad_scans:
    bad_scans_save[scankey] = scans[scankey]
    scans.pop(scankey)
    


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
            z_abs[t_name][1] += abs(t_mere["p"][2])

z_avg=[]
for t_name in z_abs:
    z_ab = z_abs[t_name]
    z_avg.append((t_name,z_ab[1]/float(z_ab[0])))

z_avg.sort(reverse = True, key = key_for_sort)
z_max_t_name = z_avg[0][0]

# having found the largest average z, we now make sure that that is positive in each
# scan. If we have to change the sign, we change the sign for all the z values in the scan.
# This performs the reflection.
scan_orientation = {}
for scankey in scans:
    scan_orientation[scankey] = "ambiguous"

def do_flip(scan):
    for t_mere in scan.values():
        t_mere["p"][2] *= -1
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"][2] *= -1
# note: I won't make reference to the positive or negative values for z_sign further

o_first = True
for z_val in z_avg: 
    t_for_reflect = z_val[0]
    o_sign = 0 
    if o_first:
        o_sign = 1
    if not o_first:
        for scankey, v in scan_orientation.items():
            if scan_orientation == "fixed" and abs(scans[scankey][t_for_reflect]["p"][2]) > 1.5:
                o_sign = sign_of(scans[scankey][t_for_reflect]["p"][2])
                break
            # we want to flip subsequent scans so that t_for_reflect ends up in the
            # same hemisphere
        if o_sign == 0:
            continue

    for scankey, scan in scans.items():
        if scan_orientation[scankey] == "ambiguous":
            if scan.get(t_for_reflect) != None:
                t_mere = scan[t_for_reflect]
                if t_mere["z_sign"] != "ambiguous" and abs(t_mere["p"][2])> 1.5:
                    if sign_of(t_mere["p"][2]) != o_sign:
                        do_flip(scan)
                    scan_orientation[scankey] = "fixed"
                    o_first = False
    all_done = True
    for scankey in scans:
        if scan_orientation[scankey] == "ambiguous":
            all_done = False
            break
    if all_done:
        break

dropcount=0
for scankey, v in scan_orientation.items():
    if v == "ambiguous":
        dropcount += 1
        scans.pop(scankey)
print(f"unable to resolve orientation for {dropcount} scans")


# bad_scans = []
# for scankey, scan in scans.items():
#     t_mere = scan.get(z_max_t_name)
#     if t_mere == None:
#         # hm.  This is a problem.  Of course this shouldn't happen
#         # but if it does we can fix it by reference to some other telomere
#         print("found a bad scan")
#         bad_scans.append(scankey)
#     elif t_mere["z_sign"] == "ambiguous" and sign_of(t_mere["p"][2]) != sign_of(t_mere["loc_alt"][2]):
#         print("found another bad scan")
#         bad_scans.append(scankey)
#     elif not sign_of(t_mere["p"][2]): # if the telomere z is negative, we have to reflect it
#         print("reflecting")
#         for t_mere in scan.values():
#             t_mere["p"][2] *= -1
#             if t_mere["z_sign"] == "ambiguous": # could happen, if a rotation put both possible locations here
#                 t_mere["loc_alt"][2] *= -1

# print(f"throwing away {len(bad_scans)} more scans")  # maybe we'll figure something better.
# for scankey in bad_scans:
#     scans.pop(scankey)


centers={}
point_versions = {}
def set_centers():
    for scan in scans.values():
        for t_name, t_mere in scan.items():
            if t_mere["z_sign"]!= "ambiguous":
                if point_versions.get(t_name) == None:
                    point_versions[t_name]=[]
                point_versions[t_name].append(t_mere["p"])

    for t_name, pvs in point_versions.items():
        centers[t_name] = utils.center_point(pvs,5) # 5 being the radius

# Now, best fit for any points left ambiguous
for scan in scans.values():
    for t_mere in scan.values():
        if t_mere["z_sign"] == "ambiguous": # find the one that's closest
            if centers.get(t_name) == None: 
                continue # nothing we can do, as they are all ambiguous
            center = centers[t_name]
            dist1 = (center[0] - t_mere["p"][0]) ** 2 + (center[1] - t_mere["p"][1]) ** 2 \
                  + (center[2] - t_mere["p"][2]) ** 2
            dist2 = (center[0] - t_mere["loc_alt"][0]) ** 2 + (center[1] - t_mere["loc_alt"][1]) ** 2 \
                  + (center[2] - t_mere["loc_alt"][2]) ** 2
            if dist1 > dist2:
                t_mere["p"] = t_mere["loc_alt"]
            t_mere["z_sign"] = "fixed"


set_centers() # one more time, with added in points
# ok, now to find the standard deviations

test_pole = telomere_names[3]
print(point_versions[test_pole])


std_deviations = {}
for t_name, values in point_versions.items():
    std_deviations[t_name] = utils.std_dev_points(point_versions[t_name], centers[t_name])

print("standard deviations")
print(std_deviations)

raise Exception("quit here")

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

    


              

            
