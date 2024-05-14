from random import random
import math


def matrix_rotate(m,t): # matrix multiplication, matrix times vector
    x = m[0][0] * t[0] + m[0][1] * t[1] + m[0][2] * t[2]
    y = m[1][0] * t[0] + m[1][1] * t[1] + m[1][2] * t[2]
    z = m[2][0] * t[0] + m[2][1] * t[1] + m[2][2] * t[2]
    return [x,y,z]

def matrix_multiply(m1, m2):
    output = [[0,0,0],[0,0,0],[0,0,0]]
    for i in range(3):
        for j in range(3):
            for k in range(3):
                output[i][j] += m1[i][k] * m2[k][j]
    return output

def random_point():
    r1 = random()*math.pi
    p = (math.sin(r1)*5.0,math.cos(r1)*5.0,0)
    r2 = random()*math.pi*2    
    y_rotate = [[math.cos(r2),0, math.sin(r2)],[0,1,0],[-math.sin(r2),0,math.cos(r2)] ]
    return matrix_rotate(y_rotate,p)

def rotate_to_x_y(y1p):
    if y1p[2] == 0: # find rotation about y axis that moves the point to the x/y plane with x >= 0
        if y1p[0] >= 0:
            return [[1,0,0],[0,1,0],[0,0,1]]
        else:
            return [[-1,0,0],[0,1,0],[0,0,-1]]
    else:
        hyp = math.sqrt(y1p[2]**2 + y1p[0]**2)
        sin_theta = y1p[2] / hyp
        cos_theta = y1p[0] / hyp
        return [[cos_theta, 0, sin_theta], [0,1,0],[-sin_theta, 0, cos_theta]]
    
def rotate_to_y_z(y1p): # find rotation about z axis that moves point to y/z plane with y >= 0
    if y1p[0] == 0:
        if y1p[1] >= 0:
            return [[1,0,0],[0,1,0],[0,0,1]]
        else:
            return [[1,0,0],[0,-1,0],[0,0,1]]
    else:
        hyp = math.sqrt(y1p[0]**2 + y1p[1]**2)
        sin_theta = y1p[0] / hyp
        cos_theta = y1p[1] / hyp
        return [[cos_theta, -sin_theta, 0], [sin_theta,cos_theta,0],[0, 0, 1]]
    
def center_point(points,radius):
    sum = [0,0,0]
    for point in points:
        for i in range(3):
            sum[i] += point[i]
    for i in range(3):
        sum[i] = sum[i]/len(points)
    length = math.sqrt(sum[0]**2 + sum[1]**2 + sum[2]**2)
    for i in range(3):
        sum[i] *= radius/length
    return sum  # keep it on the surface of the sphere
    
def std_dev_points(points,center):
    sum = 0
    for point in points:
        d = 0
        for i in range(3):
            d += (center[i] - point[i])**2
        sum += math.sqrt(d)
    sum = sum/len(points)
    return math.sqrt(sum)

def random_rotation():
    r1 = random()*math.pi
    z_rotate = [[math.cos(r1),-math.sin(r1), 0],[math.sin(r1),math.cos(r1), 0], [0,0,1]]
    r2 = random()*math.pi*2
    y_rotate = [[math.cos(r2),0, math.sin(r2)],[0,1,0],[-math.sin(r2),0,math.cos(r2)] ]
    return matrix_multiply(y_rotate,z_rotate)


    





