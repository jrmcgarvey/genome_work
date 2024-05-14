# we generate a collection of 6 points, randomly
# distributed in a sphere.  We then randomly rotate
# them.  We then rotate them back.  We then
# check the standard deviation for each

import utils

# m1 = [[1,2,3],[4,5,6],[7,8,9]]
# print(m1)
# m2 = [[1,-1,-2],[-1,1,3],[2,1,-1]]
# print(m2)
# print(utils.matrix_multiply(m1,m2))
# raise Exception("stop here")
    

points=[]
for i in range(6):
    points.append(utils.random_point())

rotations=[]
for i in range(20):
    r = utils.random_rotation()
    # print(r)
    r_points=[]
#    i = 0
    for point in points:
        r_points.append(utils.matrix_rotate(r,point))
        # if i==0:
        #     print(f"first point {utils.matrix_rotate(r,point)} was {point}")
        # i += 1
        # raise Exception("stop here")
    rotations.append(r_points)

restored = []

for r in rotations:
    r1 = utils.rotate_to_x_y(r[0])
    new_r = []
#    first_rotate = True
    for point in r:
        new_r.append(utils.matrix_rotate(r1,point))
        # if first_rotate:
        #     print(f"new z, should be 0: {utils.matrix_rotate(r1,point)[2]}")
        #     first_rotate=False      
    r2 = utils.rotate_to_y_z(new_r[0])
    new_r2 = []
#    first_rotate = True
    for point in new_r:
        new_r2.append(utils.matrix_rotate(r2,point))
        # if first_rotate:
        #     p = utils.matrix_rotate(r2,point)
        #     print(f"new x and z, should both be 0: {p[0]} {p[2]}")
        #     first_rotate = False
    r3 = utils.rotate_to_x_y(new_r2[1])
    new_r3 = []
 #   rotate_count=0
    for point in new_r2:
        new_r3.append(utils.matrix_rotate(r3,point))
        # rotate_count += 1
        # if rotate_count == 2:
        #     p = utils.matrix_rotate(r3,point)
        #     print(f"after 3rd rotate: {new_r3[0]} {new_r3[1]}")
            # raise Exception("stop here")
    restored.append(new_r3)
    # print(points)
    # print(restored[0])
    # raise Exception("stop here")


# print(restored[0][3])
# print(restored[9][3])
# raise Exception("quit here")
# restored.append(points) # add original

# print("line 40")
# print(restored[0])

# print(restored[5][3])
# print(restored[16][3])
# raise Exception("stop here")

point_versions = []
for i in range(6):
    point_versions.append([])

for r in restored:
    i = 0
    for point in r:
        point_versions[i].append(point)
        i += 1

# print(f"point versions {point_versions[3][5]} {point_versions[3][16]}")
# raise Exception("stop here")

centers = []
std_devs = []

# for point in point_versions[5]:
#     print(p[0])
# raise Exception("stop here")

# print("line 54")
# print(point_versions[0])
i = 0
for pv in point_versions:
    center = utils.center_point(pv,5)
    # print(center)
    # print(pv[0])
    # raise Exception("stop here")
    centers.append(center)
    std_devs.append(utils.std_dev_points(pv,center))

print(std_devs)  # should be an array of small values







