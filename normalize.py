import sys
import csv
import math
from copy import deepcopy
import utils

# the following settings might be adjusted according to the reliability of the data
pts_to_agree = 1
experimental_error = 0.1
# end settings

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

count = 0 
ambiguous = 0
def sign_of(n):
    if n>=0:
        return 1
    return -1
first_eliminate = True
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
        if first_eliminate:
            t_mere["z_sign"] = "positive"
        if t_mere["z_sign"] == "ambiguous":
            continue

    for t_name in telomere_names:
        distances = []
        if (t_name != refpoint):
            for scankey, scan in scans.items():
                p_ref = scan[refpoint]["p"]
                p_c = scan[t_name]["p"]
                distance1 = (p_ref[0] - p_c[0])**2 + (p_ref[1] - p_c[1])**2  + (p_ref[2] - p_c[2])**2
                distance2 = (p_ref[0] - p_c[0])**2 + (p_ref[1] - p_c[1])**2  + (p_ref[2] + p_c[2])**2
                distance3 = distance2 - distance1
                distances.append((scankey, distance3, distance1, distance2))
            distances.sort(key = key_for_sort)
            compare_distance = distances[pts_to_agree-1]
            for distance in distances[pts_to_agree-1:]:
                found_error = False
                reason = "first"
                use_long_distance = False
                if distance1 < 2 * experimental_error:
                    continue # can't tell in this case
                if distance[2] < compare_distance[2] - experimental_error:
                    use_long_distance = True
                use_short_distance = False
                if distance[3] > compare_distance[2] + experimental_error:
                    use_short_distance = True
                if (not use_long_distance) and (not use_short_distance):
                    continue # not diagnostic
                scankey = distance[0]
                scan = scans[scankey]
                t_mere = scan[t_name]
                t_mere_ref = scan[refpoint]
                if use_long_distance and use_short_distance: # should never happen
                    found_error = True
                elif use_long_distance:
                    if t_mere_ref["z_sign"] == "positive":
                        t_mere["z_sign"] = "negative"
                    else:
                        t_mere["z_sign"] = "positive"
                else:
                    t_mere["z_sign"] = t_mere_ref["z_sign"]
                if (not found_error) and includez:
                    sign1 = sign_of(t_mere["real_z"])
                    sign2 = -1
                    if t_mere["z_sign"] == "positive":
                        sign2 = 1
                    sign3 = sign_of(t_mere_ref["real_z"])
                    sign4 = -1
                    if t_mere_ref["z_sign"] == "positive":
                        sign4 = 1
                    if sign1 * sign2 * sign3 * sign4 < 0:
                        found_error = True
                        reason = "second"
                if found_error and includez:
                    print("starting error report")
                    print(f"distance {distance}")
                    print(f"compare_distance {compare_distance}")
                    print(f"t_name: {t_name} refpoint: {refpoint} reason: {reason}")
                    print(f"t_mere: {t_mere}")
                    print(f"t_mere_ref: {t_mere_ref}")
                    p1 = t_mere["p"]
                    p2 = t_mere_ref["p"]
                    s1 = sign_of(t_mere["real_z"])
                    s2 = sign_of(t_mere_ref["real_z"])
                    real_distance = (p1[0] - p2[0])**2 + (p1[1]-p2[1])**2 + (p1[2]*s1 - p2[2]*s2)**2
                    print(f"real distance {real_distance}")
                    raise Exception("Something went wrong!")
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
    if ambiguous == 0:
        break

for scan in scans.values():
    for t_mere in scan.values():
        if t_mere["z_sign"] == "negative":
            t_mere["p"][2] *= -1
        elif t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] = deepcopy(t_mere["p"])
            t_mere["loc_alt"][2] *= -1

bad_scans=[]
north_pole = telomere_names[0]
xy_plane = telomere_names[1]
for scankey, scan in scans.items():
    if scan.get(north_pole) == None:
        bad_scans.append(scankey)
        continue
    if scan.get(xy_plane) == None:
        bad_scans.append(scankey)
        continue
    if scan[xy_plane]["z_sign"] == "ambiguous":
        bad_scans.append(scankey)
        continue
    if scan[north_pole]["z_sign"] == "ambiguous": # doesn't happen, but we might change the logic
        bad_scans.append(scankey)
        continue # we might be able to put these back in, by rotating w/r to different points
    # but for now we will just toss these.
    matrix1 = utils.rotate_to_x_y(scan[north_pole]["p"])
    p_s = "p"
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix1,t_mere["p"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix1,t_mere["loc_alt"])
    matrix2 = utils.rotate_to_y_z(scan[north_pole]["p"])
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix2,t_mere["p"])

        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix2,t_mere["loc_alt"])
    matrix3 = utils.rotate_to_x_y(scan[xy_plane]["p"])
    for t_mere in scan.values():
        t_mere["p"] = utils.matrix_rotate(matrix3,t_mere["p"])
        if t_mere["z_sign"] == "ambiguous":
            t_mere["loc_alt"] == utils.matrix_rotate(matrix3, t_mere["loc_alt"])
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
            if scan_orientation == "fixed" and abs(scans[scankey][t_for_reflect]["p"][2]) > 1:
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
                if t_mere["z_sign"] != "ambiguous" and abs(t_mere["p"][2])> 1:
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

    


              

            
