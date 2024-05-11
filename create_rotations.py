import math
import random
import csv
import sys

def matrix_rotate(m,t):
    x = m[0][0] * t[0] + m[0][1] * t[1] + m[0][2] * t[2]
    y = m[1][0] * t[0] + m[1][1] * t[1] + m[1][2] * t[2]
    z = m[2][0] * t[0] + m[2][1] * t[1] + m[2][2] * t[2]
    return (x,y,z)

base = {"Y1p": (0, 5.0, 0)}
rand = random.random()*math.pi

base["Xp"] = (math.sin(rand)*4.99,math.cos(rand)*4.99,0)
for t_name in ["Xp", "Y1q", "Y2p", "Y2q"]:
    r1 = random.random()*math.pi*2
    p = (math.sin(r1)*4.99,math.cos(r1)*4.99,0)
    r2 = random.random()*math.pi*2
    y_rotate = [[math.cos(r2),0, math.sin(r2)],[0,1,0],[-math.sin(r2),0,math.cos(r2)] ]
    base[t_name] = matrix_rotate(y_rotate,p)

scans = {}
for i in range(20):
    scankey = "test" + str(i)
    scans[scankey] = {}
    scan = scans[scankey]
    r1 = random.random()*math.pi*2
    z_rotate = [[math.cos(r1),-math.sin(r1), 0],[math.sin(r1),math.cos(r1), 0], [0,0,1]]
    r2 = random.random()*math.pi*2
    y_rotate = [[math.cos(r2),0, math.sin(r2)],[0,1,0],[-math.sin(r2),0,math.cos(r2)] ]
    for t_name in base:
        scan[t_name] = matrix_rotate(y_rotate, matrix_rotate(z_rotate, base[t_name]))

f = open(sys.argv[1]+".start", "w")
for t_name in base:
    t_mere = base[t_name]
    f.write(f"{t_name} x: {t_mere[0]:.4f} y: {t_mere[1]:.4f} z: {t_mere[2]:.4f}\n")
    f.close

with open(sys.argv[1]+".csv","w") as csvfile:
    twriter = csv.writer(csvfile)
    row = [""]
    for t_name in base:
        row.extend( [t_name + "_x",t_name + "_y","","",""])
    twriter.writerow(row)
    for scankey in scans:
        scan = scans[scankey]
        row = [scankey]
        for t_name in scan:
            t_mere = scan[t_name]
            row.extend([round(t_mere[0],4), round(t_mere[1],4),"","",""])
        twriter.writerow(row)





