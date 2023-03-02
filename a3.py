'''
Ethan Hebert
10381833
2-16-23
Assignment 3
This program renders 4 objects(a pyramid, 2 cubes, and a cylinder) in 3D space with back face culling, polygon filling, 
and Z-buffering, and allows for the translation, in-place uniform scaling, and in-place rotation of each.
The first 3 render modes (1-3) switch between a wireframe, a polygon fill, or both for all objects in the scene.
The next 3 render modes (4-6) switch between flat shading, Gouraud shading, and Phong shading of the cylidner only.
Keys 7-9 are used to control the RGB values of the cylinder to view shading in different colors.
'''

# Import libraries
import math
import copy
from tkinter import *

# Constants
CANVASWIDTH = 400
CANVASHEIGHT = 400
D = 500
SELECTED_EDGE_COLOR = "#00fff7"
CYLINDER_SCENE_NO = 3

# Lighting constants
IA = 0.25 # intensity of the ambient light in the scene
IP = 1 - IA # intensity of the point light source in the scene
L = [1,1,-1] # lighting vector, 45 degree angle, light behind viewer's right shoulder
V = [0,0,-1] # view vector, points towards viewer/center of projection
KD = 0.75 # constant of diffuse reflectivity
KS = 1 - KD # constant of specular reflectivity
SPECINDEX = 2 # specular index, the spread of relfected light in the mirrored direction
AMBIENT = IA*KD # ambient diffuse component of illumination model - constant for every polygon

# ***************************** Initialize Pyramid Object ***************************
# Definition  of the five underlying points
p_apex = [150,-100,100]
p_base1 = [200,-200,50]
p_base2 = [200,-200,150]
p_base3 = [100,-200,150]
p_base4 = [100,-200,50]

# Definition of the five polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
p_frontpoly = [p_apex,p_base1,p_base4]
p_rightpoly = [p_apex,p_base2,p_base1]
p_backpoly = [p_apex,p_base3,p_base2]
p_leftpoly = [p_apex,p_base4,p_base3]
p_bottompoly = [p_base1,p_base2,p_base3,p_base4]

# Definition of the object
Pyramid = [p_bottompoly, p_frontpoly, p_rightpoly, p_backpoly, p_leftpoly]
# Defintion of the color associated with each face of the pyramid
PyramidColor = ["black", "red", "green", "blue", "yellow"]

# Definition of the Pyramid's underlying point cloud.  No structure, just the points.
PyramidPointCloud = [p_apex, p_base1, p_base2, p_base3, p_base4]
DefaultPyramidPointCloud = copy.deepcopy(PyramidPointCloud)
#************************************************************************************

# ***************************** Initialize Cube 1 ***************************
# Definition  of the eight underlying points
c1_top1 = [-100,150,50]
c1_top2 = [-100,150,100]
c1_top3 = [-150,150,100]
c1_top4 = [-150,150,50]
c1_base1 = [-100,100,50]
c1_base2 = [-100,100,100]
c1_base3 = [-150,100,100]
c1_base4 = [-150,100,50]

# Definition of the six polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
c1_frontpoly = [c1_top1,c1_base1,c1_base4,c1_top4]
c1_rightpoly = [c1_top2,c1_base2,c1_base1,c1_top1]
c1_backpoly = [c1_top3,c1_base3,c1_base2,c1_top2]
c1_leftpoly = [c1_top4,c1_base4,c1_base3,c1_top3]
c1_bottompoly = [c1_base1,c1_base2,c1_base3,c1_base4]
c1_toppoly = [c1_top4,c1_top3,c1_top2,c1_top1]

# Definition of the object
Cube1 = [c1_frontpoly,c1_backpoly,c1_rightpoly,c1_leftpoly,c1_bottompoly,c1_toppoly]
# Defintion of the color associated with each face of the cube
Cube1Color = ["white", "#cccccc", "#999999", "#666666","#333333", "black"]

# Definition of Cube1's underlying point cloud.  No structure, just the points.
Cube1PointCloud = [c1_top1,c1_top2,c1_top3,c1_top4,c1_base1,c1_base2,c1_base3,c1_base4]
DefaultCube1PointCloud = copy.deepcopy(Cube1PointCloud)
#************************************************************************************

# ***************************** Initialize Cube 2 ***************************
# Definition  of the eight underlying points
c2_top1 = [150,150,50]
c2_top2 = [150,150,100]
c2_top3 = [100,150,100]
c2_top4 = [100,150,50]
c2_base1 = [150,100,50]
c2_base2 = [150,100,100]
c2_base3 = [100,100,100]
c2_base4 = [100,100,50]

# Definition of the six polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
c2_frontpoly = [c2_top1,c2_base1,c2_base4,c2_top4]
c2_rightpoly = [c2_top2,c2_base2,c2_base1,c2_top1]
c2_backpoly = [c2_top3,c2_base3,c2_base2,c2_top2]
c2_leftpoly = [c2_top4,c2_base4,c2_base3,c2_top3]
c2_bottompoly = [c2_base1,c2_base2,c2_base3,c2_base4]
c2_toppoly = [c2_top4,c2_top3,c2_top2,c2_top1]

# Definition of the object
Cube2 = [c2_frontpoly,c2_backpoly,c2_rightpoly,c2_leftpoly,c2_bottompoly,c2_toppoly]
# Defintion of the color associated with each face of the cube
Cube2Color = ["white", "#cccccc", "#999999", "#666666","#333333", "black"]

# Definition of Cube2's underlying point cloud.  No structure, just the points.
Cube2PointCloud = [c2_top1,c2_top2,c2_top3,c2_top4,c2_base1,c2_base2,c2_base3,c2_base4]
DefaultCube2PointCloud = copy.deepcopy(Cube2PointCloud)
#************************************************************************************

# ***************************** Initialize Cylinder Object ***************************
# Definition of the 16 underlying points
cyl_front1 = [-25,60.35535,50]
cyl_front2 = [25,60.35535,50]
cyl_front3 = [60.35535,25,50]
cyl_front4 = [60.35535,-25,50]
cyl_front5 = [25,-60.35535,50]
cyl_front6 = [-25,-60.35535,50]
cyl_front7 = [-60.35535,-25,50]
cyl_front8 = [-60.35535,25,50]
cyl_back1 = [-25,60.35535,250]
cyl_back2 = [25,60.35535,250]
cyl_back3 = [60.35535,25,250]
cyl_back4 = [60.35535,-25,250]
cyl_back5 = [25,-60.35535,250]
cyl_back6 = [-25,-60.35535,250]
cyl_back7 = [-60.35535,-25,250]
cyl_back8 = [-60.35535,25,250]

# Definition of the ten polygon faces using the meaningful point names
# Polys are defined in clockwise order when viewed from the outside
cyl_northPoly = [cyl_front1, cyl_back1, cyl_back2, cyl_front2]
cyl_northEastPoly = [cyl_front2, cyl_back2, cyl_back3, cyl_front3]
cyl_eastPoly = [cyl_front3, cyl_back3, cyl_back4, cyl_front4]
cyl_southEastPoly = [cyl_front4, cyl_back4, cyl_back5, cyl_front5]
cyl_southPoly = [cyl_front5, cyl_back5, cyl_back6, cyl_front6]
cyl_southWestPoly = [cyl_front6, cyl_back6, cyl_back7, cyl_front7]
cyl_westPoly = [cyl_front7, cyl_back7, cyl_back8, cyl_front8]
cyl_northWestPoly = [cyl_front8, cyl_back8, cyl_back1, cyl_front1]
cyl_frontPoly = [cyl_front1, cyl_front2, cyl_front3, cyl_front4, cyl_front5, cyl_front6, cyl_front7, cyl_front8]
cyl_backPoly = [cyl_back1, cyl_back8, cyl_back7, cyl_back6, cyl_back5, cyl_back4, cyl_back3, cyl_back2]

# Definition of the cylinder object
Cylinder = [cyl_northPoly, cyl_northEastPoly, cyl_eastPoly, cyl_southEastPoly, cyl_southPoly, 
            cyl_southWestPoly, cyl_westPoly, cyl_northWestPoly, cyl_frontPoly, cyl_backPoly]
# Defintion of the color associated with each face of the cylinder
CylinderColor = ["#d3fcd5","#44344f","#c2f970","#564d80","#d3fcd5","#44344f","#c2f970","#564d80","#98a6d4","#98a6d4"]

# Definition of Cylinder's underlying point cloud.  No structure, just the points.
CylinderPointCloud = [cyl_front1, cyl_front2, cyl_front3, cyl_front4, cyl_front5, cyl_front6, cyl_front7, cyl_front8,
                      cyl_back1, cyl_back2, cyl_back3, cyl_back4, cyl_back5, cyl_back6, cyl_back7, cyl_back8]
DefaultCylinderPointCloud = copy.deepcopy(CylinderPointCloud)
#************************************************************************************

# Arrays to hold all objects in the scene to allow switching selected object.
# Variable to store the index of the currently selected object as well as render mode.
scene = [Cube1,Pyramid,Cube2,Cylinder]
scenePointClouds = [Cube1PointCloud,PyramidPointCloud,Cube2PointCloud,CylinderPointCloud]
sceneDefaultPointClouds = [DefaultCube1PointCloud,DefaultPyramidPointCloud,DefaultCube2PointCloud,DefaultCylinderPointCloud]
sceneColor = [Cube1Color,PyramidColor,Cube2Color,CylinderColor]
currObject = 1
renderMode = 1 # 1 = wireframe (edges) | 2 = polyfill + wireframe | 3 = polyfill
               # 4 = flat shading | 5 = Gouraud shading | 6 = Phong shading (4-6 are cylinder only, no wireframe)
rgb = [False,True,False] # bools to allow switching color of cylinder using keys 7-9 in renderModes 4-6

# This function resets an object to its original size and location in 3D space.
def resetObject(objectIndex):
    pointCloud = scenePointClouds[objectIndex]
    defaultPointCloud = sceneDefaultPointClouds[objectIndex]
    for i in range(len(pointCloud)):
        for j in range(3):
            pointCloud[i][j] = defaultPointCloud[i][j]

# This function returns the reference point (visual center) of an object.
# Finds the max and min x, y, and z values and returs the average of each.
def referencePoint(object):
    max = [object[0][0], object[0][1], object[0][2]] # first point as default values
    min = [object[0][0], object[0][1], object[0][2]] # to compare to

    for i in range(len(object)): # find max and min
        for j in range(3):
            if (object[i][j] > max[j]):
                max[j] = object[i][j]
            if (object[i][j] < min[j]):
                min[j] = object[i][j]

    centerx = (min[0]+max[0])/2 # return reference point
    centery = (min[1]+max[1])/2
    centerz = (min[2]+max[2])/2
    refPoint = [centerx, centery, centerz]
    return refPoint

# This function translates an object by some displacement.  The displacement is a 3D
# vector so the amount of displacement in each dimension can vary.
def translate(object,displacement):
    for i in range(len(object)):
        for j in range(3):
            object[i][j] += displacement[j] # translate each point by its displacement
    
# This function performs a simple uniform scale of an object assuming the object is
# centered at the origin.  The scalefactor is a scalar.
def scale(object,scalefactor):
    refPoint = referencePoint(object) # translate reference point to origin
    translate(object, [-refPoint[0], -refPoint[1], -refPoint[2]])

    for i in range(len(object)): # scale
        for j in range(3):
            object[i][j] *= scalefactor

    translate(object,refPoint) #translate reference point back to original position

# This function converts from degrees to radians. 
def degToRad(deg):
    rad = deg*math.pi/180
    return rad

# This function performs a rotation of an object about the Z axis (from +X to +Y)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CCW
# in a LHS when viewed from -Z [the location of the viewer in the standard postion]
def rotateZ(object,degrees):
    refPoint = referencePoint(object) # translate reference point to origin
    translate(object, [-refPoint[0], -refPoint[1], -refPoint[2]])

    phi = degToRad(degrees) # rotate
    for i in range(len(object)): # rotate, z value doesn't change
        x = object[i][0]
        y = object[i][1]
        xrz = x*math.cos(phi) - y*math.sin(phi)
        yrz = x*math.sin(phi) + y*math.cos(phi)
        object[i][0] = xrz
        object[i][1] = yrz

    translate(object,refPoint) #translate reference point back to original position
    
# This function performs a rotation of an object about the Y axis (from +Z to +X)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CW
# in a LHS when viewed from +Y looking toward the origin.
def rotateY(object,degrees):
    refPoint = referencePoint(object) # translate reference point to origin
    translate(object, [-refPoint[0], -refPoint[1], -refPoint[2]])

    phi = degToRad(degrees) # rotate
    for i in range(len(object)): # rotate, y value doesn't change
        x = object[i][0]
        z = object[i][2]
        xry = x*math.cos(phi) + z*math.sin(phi)
        zry = -x*math.sin(phi) + z*math.cos(phi)
        object[i][0] = xry
        object[i][2] = zry

    translate(object,refPoint) #translate reference point back to original position

# This function performs a rotation of an object about the X axis (from +Y to +Z)
# by 'degrees', assuming the object is centered at the origin.  The rotation is CW
# in a LHS when viewed from +X looking toward the origin.
def rotateX(object,degrees):
    refPoint = referencePoint(object) # translate reference point to origin
    translate(object, [-refPoint[0], -refPoint[1], -refPoint[2]])

    phi = degToRad(degrees) # rotate
    for i in range(len(object)): # rotate, x value doesn't change
        y = object[i][1]
        z = object[i][2]
        yrx = y*math.cos(phi) - z*math.sin(phi)
        zrx = y*math.sin(phi) + z*math.cos(phi)
        object[i][1] = yrx
        object[i][2] = zrx

    translate(object,refPoint) #translate reference point back to original position

# This function will draw the scene by repeatedly calling drawObject on each object in the scene.
# If the object is the currently selected object, "True" is passed down the functions to eventually
# draw the edges a color other than black and bold to show selection. 
# Each object's polygon face colors are also passed down. renderMode is passed down as well.
# Here the z buffer is initlaized to all values of D.
# renderModes 4-6 call their respective shading functions with the cylinder object
def drawScene(scene,currObject,renderMode):
    # zBuffer array with the z depth of every pixel on the screen.
    # Default value is D, array size is size of canvas.
    zBuffer = [[D for y in range(CANVASHEIGHT)] for x in range(CANVASWIDTH)]

    # shading of cylinder only
    # flat shading (renderMode 4)
    if (renderMode == 4):
        flatShading(scene[CYLINDER_SCENE_NO],zBuffer)
    # gourard shading (renderMode 5)
    elif (renderMode == 5):
        gouraudShading(scene[CYLINDER_SCENE_NO],zBuffer)
    # phong shading (renderMode 6)
    elif (renderMode == 6):
        phongShading(scene[CYLINDER_SCENE_NO],zBuffer)

    # polygon fill of all objects with z buffer
    else:
        for i in range(len(scene)):
            if (i == currObject):
                drawObject(scene[i],sceneColor[i],True,renderMode,zBuffer) # draw in blue edges
            else:
                drawObject(scene[i],sceneColor[i],False,renderMode,zBuffer) # draw in black edges

# This function will draw an object by repeatedly callying drawPoly on each polygon in the object.
def drawObject(object,objectColor,isCurrObject,renderMode,zBuffer):
    for i in range(len(object)):
        drawPoly(object[i],objectColor[i],isCurrObject,renderMode,zBuffer)

# This function will draw a polygon by repeatedly callying drawLine on each pair of points
# making up the object. It will only draw the polygon if it is visible after back face culling.
# Each polygon is drawn its predefined color from each object's color array.
# This function uses the renderMode to draw the edges, polygon fill, or both.
def drawPoly(poly,polyColor,isCurrObject,renderMode,zBuffer):
    if (backFaceCulling(poly)): # only draw the polygon if it is visible after back face culling
        # normal polygon fill (renderMode 2 and 3)
        if (renderMode == 2 or renderMode == 3):
            polyFill(poly,polyColor,zBuffer)
        # draw the edges (renderMode 1 and 2)
        if (renderMode == 1 or renderMode == 2):
            for i in range(len(poly)):
                drawLine(poly[i], poly[(i+1)%len(poly)],isCurrObject)

# This function draws the line using the built-in create_line method. currObject is drawn in bold in a different color than black.
def drawLine(start,end,isCurrObject):
    startdisplay = convertToDisplayCoordinates(project(start)) #convert 3D -> 2D -> display
    enddisplay = convertToDisplayCoordinates(project(end))
    if (isCurrObject): # if this is currObject, draw the lines red and bold
        w.create_line(startdisplay[0], startdisplay[1], enddisplay[0], enddisplay[1], width=3, fill=SELECTED_EDGE_COLOR)
    else:
        w.create_line(startdisplay[0], startdisplay[1], enddisplay[0], enddisplay[1])

# This function converts from 3D to 2D (+ depth) using the perspective projection technique.
def project(point):
    x = point[0]
    y = point[1]
    z = point[2]
    ps = [D*x/(D+z), D*y/(D+z), D*z/(D+z)]
    return ps

# This function converts a 2D point to display coordinates in the tk system with an included z depth value.
def convertToDisplayCoordinates(point):
    xps = point[0]
    yps = point[1]
    zps = point[2]
    displayXYZ = [(CANVASWIDTH/2)+xps, (CANVASHEIGHT/2)-yps, zps]
    return displayXYZ
    
# This function calculates the surface normal of a polygon
def surfaceNormal(poly):
    # create 2 "vectors" p and q which go from P0 to P1 and P0 to P2 (points) on the polygon.
    p = [poly[1][0]-poly[0][0], poly[1][1]-poly[0][1], poly[1][2]-poly[0][2]]
    q = [poly[2][0]-poly[0][0], poly[2][1]-poly[0][1], poly[2][2]-poly[0][2]]

    # calculate the surface normal n of the polygon using p and q (cross product).
    n = [p[1]*q[2]-p[2]*q[1], p[2]*q[0]-p[0]*q[2], p[0]*q[1]-p[1]*q[0]]
    return n

# Convert an N dimensional vector into a unit vector (normalize the vector)
def normalize(vector):
    sumOfSquares = 0
    for i in range(len(vector)):
        sumOfSquares += vector[i]**2
    magnitude = math.sqrt(sumOfSquares)
    vect = []
    for i in range(len(vector)):
        vect.append(vector[i]/magnitude)
    return vect

# This function evaluates if a polygon faces the viewport or not to determine if it should be drawn or not.
def backFaceCulling(poly):
    # get the normalized surface normal of the polygon
    nNorm = normalize(surfaceNormal(poly))

    # compute the plane offset (determine the position of the plane containing the polygon - plane equation)
    d = nNorm[0]*poly[0][0] + nNorm[1]*poly[0][1] + nNorm[2]*poly[0][2]

    # determine visibility of the polygon from viewpoint v = <0,0,-D> using dot product
    v = [0,0,-D]
    visTest = nNorm[0]*v[0] + nNorm[1]*v[1] + nNorm[2]*v[2] - d
    if (visTest > 0): # viewpoint is above the face of the polygon, visible
        return True
    else: # viewpoint is inside of or behind the face of the polygon, invisible
        return False

# This function fills in a polygon face a given color to make it appear solid.
# Z buffer is also implemented here.
def polyFill(poly,color,zBuffer):
    # convert polygon vertices to projected display coordinates, rounded reals (not ints)
    # don't round the z value to compare depths with total accuracy
    polyDisplay = []
    for i in range(len(poly)):
        polyDisplay.append(convertToDisplayCoordinates(project(poly[i])))
        polyDisplay[i][0] = round(polyDisplay[i][0],0)
        polyDisplay[i][1] = round(polyDisplay[i][1],0)

    # pre-compute the edge constants
    edgeTable = computeEdgeTable(polyDisplay)

    # if edge table is empty, then no need to fill anything! exit function now.
    if not edgeTable:
        return

    # calculate the range of y values that the fill lines will draw between
    firstFillLine = edgeTable[0]["yStart"]
    yEndVals = []
    for i in range(len(edgeTable)):
        yEndVals.append(edgeTable[i]["yEnd"])
    lastFillLine = max(yEndVals)

    # indices for the first (i), second (j), and next (next) edges
    i = 0
    j = 1
    next = 2

    # find the xStart and zStart values on each edge to know where to draw the first line btwn
    ix = edgeTable[i]["xStart"]
    jx = edgeTable[j]["xStart"]
    iz = edgeTable[i]["zStart"]
    jz = edgeTable[j]["zStart"]

    # paint one fill line at a time (y)
    for y in range(int(firstFillLine), int(lastFillLine)+1):
        # determine the left and right edge
        if (ix < jx):
            leftx = ix
            leftz = iz
            rightx = jx
            rightz = jz
        else:
            leftx = jx
            leftz = jz
            rightx = ix
            rightz = iz

        # the initial z for the current fill line
        z = leftz
        
        # compute dZ for the fill line. can be 0 if line is 1 pixel long
        if (rightx-leftx != 0):
            dZFillLine = (rightz-leftz)/(rightx-leftx)
        else:
            dZFillLine = 0

        # paint across a fill line from left to right
        for x in range(int(leftx), int(rightx)+1):
            if (0 < x < CANVASWIDTH and 0 < y < CANVASHEIGHT): # handles index out of bounds error
                if z < zBuffer[x][y]:
                    w.create_line(x,y,x+1,y,fill=color) # sets a pixel
                    zBuffer[x][y] = z
            z += dZFillLine

        # update the x and z values of edges i and j for the next fill line (add dX and dZ)
        ix += edgeTable[i]["dX"]
        jx += edgeTable[j]["dX"]
        iz += edgeTable[i]["dZ"]
        jz += edgeTable[j]["dZ"]

        # if reached the bottom of an edge, switch out to next edge until reach bottom
        if (y >= edgeTable[i]["yEnd"] and y < lastFillLine):
            i = next
            ix = edgeTable[i]["xStart"]
            iz = edgeTable[i]["zStart"]
            next += 1
        if (y >= edgeTable[j]["yEnd"] and y < lastFillLine):
            j = next
            jx = edgeTable[j]["xStart"]
            jz = edgeTable[j]["zStart"]
            next += 1

# This function pre-computes the edge table to be used in polygon filling
def computeEdgeTable(polyDisplay):
    # create a list of all the edges on the polygon, defined from point w/
    # smaller y value to point w/ larger y value, if the y values are equal
    # then don't add the edge - horizontal line
    edges = []
    for i in range(len(polyDisplay)):
        ptA = polyDisplay[i]
        ptB = polyDisplay[(i+1)%len(polyDisplay)]
        if (ptA[1] < ptB[1]):
            edges.append([ptA, ptB])
        elif (ptA[1] > ptB[1]):
            edges.append([ptB, ptA])
    
    # sort the edges in order of increasing yStart values
    # key is a function to sort by: a lambda function to find the yStart value of each edge
    edges.sort(key=(lambda edge : edge[0][1]))

    # create the final edge table, an array of dictionaries which store values and their var names
    edgeTable = []
    for i in range(len(edges)):
        edgeTable.append({
            "edge": edges[i],
            "xStart": edges[i][0][0],
            "yStart": edges[i][0][1],
            "yEnd": edges[i][1][1],
            "dX": (edges[i][1][0]-edges[i][0][0])/(edges[i][1][1]-edges[i][0][1]),
            "zStart": edges[i][0][2],
            "dZ": (edges[i][1][2]-edges[i][0][2])/(edges[i][1][1]-edges[i][0][1])
        })

    return edgeTable

# This function runs phong illumination on a vector and returns the diffuse and specular reflectivities (ambient is constant)
def phongIllumination(n):
    # normalize lighting and view vectors
    l = normalize(L)
    v = normalize(V)
    # find diffuse reflectivity
    ndotl = n[0]*l[0] + n[1]*l[1] + n[2]*l[2]
    if ndotl < 0: ndotl = 0
    diffuse = IP*KD*ndotl
    # find specular reflectivity
    r = reflect(n,l) # return vector is normalized in "reflect"
    rdotv = r[0]*v[0] + r[1]*v[1] + r[2]*v[2]
    if rdotv < 0: rdotv = 0
    specular = IP*KS*rdotv**SPECINDEX
    # return the ambient, diffuse, and specular values in an array
    return([diffuse,specular])

# This function takes in a Cylinder and returns its polygon surface normals
def cylSurfaceNormals(Cylinder):
    # find the normalized surface normal of each polygon (n)
    n = []
    for i in range(len(Cylinder)):
        n.append(normalize(surfaceNormal(Cylinder[i])))
    return n

# This function takes in a Cylinder and its polygon surface normals and returns every vertex normal
def cylVertexNormals(Cylinder,n):
    # polygons is an array to hold each polygon as a set of vertices
    polygons = []
    for i in range(len(Cylinder)):
        # vertices is an array to hold the vertices in a single polygon
        vertices = []
        for j in range(len(Cylinder[i])):
            vertices.append(Cylinder[i][j])
        polygons.append(vertices)

    vertexNormalPolys = []
    # find the surface normals at each vertex of each polygon
    for i in range(len(polygons)):
        vertexNormals = []
        # for the octagon end caps, every vertex has the polygon's surface normal
        if (i == 8 or i == 9): 
            for j in range(8):
                vertexNormals.append(n[i])
        # for the normal square side polygons, normalize the sum of touching polygon faces to get vertex normal
        else:
            # each polygon has 4 vertices
            for j in range(4):
                vector = []
                # each vector has an x, y, and z
                for k in range(3):
                    # see if you're on the 2 left-side points or 2 right-side points to determine which side polygon to add with
                    # left
                    if (j == 0 or j == 1):
                        vector.append(n[i][k] + n[(i-1)%(len(polygons)-2)][k])
                    # right
                    else:
                        vector.append(n[i][k] + n[(i+1)%(len(polygons)-2)][k])
                # append the normal of the sum vector you found to the vertexNormals array
                vertexNormals.append(normalize(vector))

        # append all the vertex normals for this poly to vertexNormalPolys
        vertexNormalPolys.append(vertexNormals)
    
    return vertexNormalPolys

# This function renders the cylinder in the flat shading method
def flatShading(Cylinder,zBuffer):
    # find the normalized surface normal of each polygon (n)
    n = cylSurfaceNormals(Cylinder)

    # find the phong illumination and color of each polygon based on its surface normal
    colors = []
    for i in range(len(Cylinder)):
        intensity = phongIllumination(n[i])
        colors.append(triColorHexCode(AMBIENT,intensity[0],intensity[1]))

    # draw each poly with its color if it passes backface culling
    for i in range(len(Cylinder)):
        if (backFaceCulling(Cylinder[i])):
            polyFill(Cylinder[i],colors[i],zBuffer)

# This function renders the cylinder in the Gouraud shading method
def gouraudShading(Cylinder,zBuffer):
    # find the normalized surface normal of each polygon (n)
    n = cylSurfaceNormals(Cylinder)

    # find the vertex normals for the Cylinder
    vertexNormalPolys = cylVertexNormals(Cylinder,n)

    # find the phong illumination at each vertex of each polygon
    intensityPolys = []
    for i in range(len(vertexNormalPolys)):
        intensityVertices = []
        for j in range(len(vertexNormalPolys[i])):
            # append the phong illumination for each vertex
            intensityVertices.append(phongIllumination(vertexNormalPolys[i][j]))
        # append the phong illuminations for this polygon
        intensityPolys.append(intensityVertices)

    # shade in each polygon with a bi-linear interpolation of intensity values if it passes backface culling
    for i in range(len(Cylinder)):
        if (backFaceCulling(Cylinder[i])):
            polyGouraudShade(Cylinder[i],intensityPolys[i],zBuffer)

# This function fills in a polygon with varying color in Gouraud shading.
def polyGouraudShade(poly,polyIntensity,zBuffer):
    # convert polygon vertices to projected display coordinates, rounded reals (not ints)
    # don't round the z value to compare depths with total accuracy
    polyDisplay = []
    for i in range(len(poly)):
        polyDisplay.append(convertToDisplayCoordinates(project(poly[i])))
        polyDisplay[i][0] = round(polyDisplay[i][0],0)
        polyDisplay[i][1] = round(polyDisplay[i][1],0)

    # pre-compute the edge constants
    edgeTable = computeEdgeTableGouraudShade(polyDisplay,polyIntensity)

    # if edge table is empty, then no need to fill anything! exit function now.
    if not edgeTable:
        return

    # calculate the range of y values that the fill lines will draw between
    firstFillLine = edgeTable[0]["yStart"]
    yEndVals = []
    for i in range(len(edgeTable)):
        yEndVals.append(edgeTable[i]["yEnd"])
    lastFillLine = max(yEndVals)

    # indices for the first (i), second (j), and next (next) edges
    i = 0
    j = 1
    next = 2

    # find the xStart and zStart values on each edge to know where to draw the first line btwn
    # also d and s for diffuse and specular intensity values
    ix = edgeTable[i]["xStart"]
    jx = edgeTable[j]["xStart"]
    iz = edgeTable[i]["zStart"]
    jz = edgeTable[j]["zStart"]
    id = edgeTable[i]["dStart"]
    jd = edgeTable[j]["dStart"]
    iS = edgeTable[i]["sStart"]
    jS = edgeTable[j]["sStart"]

    # paint one fill line at a time (y)
    for y in range(int(firstFillLine), int(lastFillLine)+1):
        # determine the left and right edge
        if (ix < jx):
            leftx = ix
            leftz = iz
            leftd = id
            lefts = iS
            rightx = jx
            rightz = jz
            rightd = jd
            rights = jS
        else:
            leftx = jx
            leftz = jz
            leftd = jd
            lefts = jS
            rightx = ix
            rightz = iz
            rightd = id
            rights = iS

        # the initial z, d, and s for the current fill line
        z = leftz
        d = leftd
        s = lefts
        
        # compute dZ, dD, and dS for the fill line. can be 0 if line is 1 pixel long
        if (rightx-leftx != 0):
            dZFillLine = (rightz-leftz)/(rightx-leftx)
            dDFillLine = (rightd-leftd)/(rightx-leftx)
            dSFillLine = (rights-lefts)/(rightx-leftx)
        else:
            dZFillLine = 0
            dDFillLine = 0
            dSFillLine = 0

        # paint across a fill line from left to right
        for x in range(int(leftx), int(rightx)+1):
            if (0 < x < CANVASWIDTH and 0 < y < CANVASHEIGHT): # handles index out of bounds error
                if z < zBuffer[x][y]:
                    w.create_line(x,y,x+1,y,fill=triColorHexCode(AMBIENT,d,s)) # sets a pixel at its phong color
                    zBuffer[x][y] = z
            z += dZFillLine
            d += dDFillLine
            s += dSFillLine

        # update the x and z values of edges i and j for the next fill line (add dX and dZ)
        # do the same for d and s
        ix += edgeTable[i]["dX"]
        jx += edgeTable[j]["dX"]
        iz += edgeTable[i]["dZ"]
        jz += edgeTable[j]["dZ"]
        id += edgeTable[i]["dD"]
        jd += edgeTable[j]["dD"]
        iS += edgeTable[i]["dS"]
        jS += edgeTable[j]["dS"]

        # if reached the bottom of an edge, switch out to next edge until reach bottom
        if (y >= edgeTable[i]["yEnd"] and y < lastFillLine):
            i = next
            ix = edgeTable[i]["xStart"]
            iz = edgeTable[i]["zStart"]
            id = edgeTable[i]["dStart"]
            iS = edgeTable[i]["sStart"]
            next += 1
        if (y >= edgeTable[j]["yEnd"] and y < lastFillLine):
            j = next
            jx = edgeTable[j]["xStart"]
            jz = edgeTable[j]["zStart"]
            jd = edgeTable[j]["dStart"]
            jS = edgeTable[j]["sStart"]
            next += 1

# This function pre-computes the edge table to be used in Gouraud polygon shading (intensity values)
def computeEdgeTableGouraudShade(polyDisplay,polyIntensity):
    # create a list of all the edges on the polygon, defined from point w/
    # smaller y value to point w/ larger y value, if the y values are equal
    # then don't add the edge - horizontal line
    edges = []
    for i in range(len(polyDisplay)):
        ptA = polyDisplay[i]
        ptB = polyDisplay[(i+1)%len(polyDisplay)]
        iA = polyIntensity[i]
        iB = polyIntensity[(i+1)%len(polyIntensity)]
        if (ptA[1] < ptB[1]):
            edges.append([ptA, ptB, iA, iB])
        elif (ptA[1] > ptB[1]):
            edges.append([ptB, ptA, iB, iA])
    
    # sort the edges in order of increasing yStart values
    # key is a function to sort by: a lambda function to find the yStart value of each edge
    edges.sort(key=(lambda edge : edge[0][1]))

    # create the final edge table, an array of dictionaries which store values and their var names
    # d is diffuse, s is specular - intensity
    edgeTable = []
    for i in range(len(edges)):
        edgeTable.append({
            "edge": edges[i],
            "xStart": edges[i][0][0],
            "yStart": edges[i][0][1],
            "yEnd": edges[i][1][1],
            "dX": (edges[i][1][0]-edges[i][0][0])/(edges[i][1][1]-edges[i][0][1]),
            "zStart": edges[i][0][2],
            "dZ": (edges[i][1][2]-edges[i][0][2])/(edges[i][1][1]-edges[i][0][1]),
            "dStart": edges[i][2][0],
            "dD": (edges[i][3][0]-edges[i][2][0])/(edges[i][1][1]-edges[i][0][1]),
            "sStart": edges[i][2][1],
            "dS": (edges[i][3][1]-edges[i][2][1])/(edges[i][1][1]-edges[i][0][1])
        })

    return edgeTable

# This function renders the cylinder in the Phong shading method
def phongShading(Cylinder,zBuffer):
    # find the normalized surface normal of each polygon (n)
    n = cylSurfaceNormals(Cylinder)

    # find the vertex normals for the Cylinder
    vertexNormalPolys = cylVertexNormals(Cylinder,n)

    # shade in each polygon with a bi-linear interpolation of vertex normals if it passes backface culling
    for i in range(len(Cylinder)):
        if (backFaceCulling(Cylinder[i])):
            polyPhongShade(Cylinder[i],vertexNormalPolys[i],zBuffer)

# This function fills in a polygon with varying color in Phong shading.
def polyPhongShade(poly,polyVertexNormals,zBuffer):
    # convert polygon vertices to projected display coordinates, rounded reals (not ints)
    # don't round the z value to compare depths with total accuracy
    polyDisplay = []
    for i in range(len(poly)):
        polyDisplay.append(convertToDisplayCoordinates(project(poly[i])))
        polyDisplay[i][0] = round(polyDisplay[i][0],0)
        polyDisplay[i][1] = round(polyDisplay[i][1],0)

    # pre-compute the edge constants
    edgeTable = computeEdgeTablePhongShade(polyDisplay,polyVertexNormals)

    # if edge table is empty, then no need to fill anything! exit function now.
    if not edgeTable:
        return

    # calculate the range of y values that the fill lines will draw between
    firstFillLine = edgeTable[0]["yStart"]
    yEndVals = []
    for i in range(len(edgeTable)):
        yEndVals.append(edgeTable[i]["yEnd"])
    lastFillLine = max(yEndVals)

    # indices for the first (i), second (j), and next (next) edges
    i = 0
    j = 1
    next = 2

    # find the xStart and zStart values on each edge to know where to draw the first line btwn
    # also nx, ny, nz for surface normals
    ix = edgeTable[i]["xStart"]
    jx = edgeTable[j]["xStart"]
    iz = edgeTable[i]["zStart"]
    jz = edgeTable[j]["zStart"]
    inx = edgeTable[i]["nxStart"]
    jnx = edgeTable[j]["nxStart"]
    iny = edgeTable[i]["nyStart"]
    jny = edgeTable[j]["nyStart"]
    inz = edgeTable[i]["nzStart"]
    jnz = edgeTable[j]["nzStart"]

    # paint one fill line at a time (y)
    for y in range(int(firstFillLine), int(lastFillLine)+1):
        # determine the left and right edge
        if (ix < jx):
            leftx = ix
            leftz = iz
            leftnx = inx
            leftny = iny
            leftnz = inz
            rightx = jx
            rightz = jz
            rightnx = jnx
            rightny = jny
            rightnz = jnz
        else:
            leftx = jx
            leftz = jz
            leftnx = jnx
            leftny = jny
            leftnz = jnz
            rightx = ix
            rightz = iz
            rightnx = inx
            rightny = iny
            rightnz = inz

        # the initial z, nx, ny, and nz for the current fill line
        z = leftz
        nx = leftnx
        ny = leftny
        nz = leftnz
        
        # compute dZ, dNX, dNY, and dNZ for the fill line. can be 0 if line is 1 pixel long
        if (rightx-leftx != 0):
            dZFillLine = (rightz-leftz)/(rightx-leftx)
            dNXFillLine = (rightnx-leftnx)/(rightx-leftx)
            dNYFillLine = (rightny-leftny)/(rightx-leftx)
            dNZFillLine = (rightnz-leftnz)/(rightx-leftx)
        else:
            dZFillLine = 0
            dNXFillLine = 0
            dNYFillLine = 0
            dNZFillLine = 0

        # paint across a fill line from left to right
        for x in range(int(leftx), int(rightx)+1):
            if (0 < x < CANVASWIDTH and 0 < y < CANVASHEIGHT): # handles index out of bounds error
                if z < zBuffer[x][y]:
                    intensity = phongIllumination(normalize([nx,ny,nz])) # get the intensity values for this point's unique surface normal
                    w.create_line(x,y,x+1,y,fill=triColorHexCode(AMBIENT,intensity[0],intensity[1])) # sets a pixel at its phong color
                    zBuffer[x][y] = z
            z += dZFillLine
            nx += dNXFillLine
            ny += dNYFillLine
            nz += dNZFillLine

        # update the x and z values of edges i and j for the next fill line (add dX and dZ)
        # do the same for nx, ny, and nz
        ix += edgeTable[i]["dX"]
        jx += edgeTable[j]["dX"]
        iz += edgeTable[i]["dZ"]
        jz += edgeTable[j]["dZ"]
        inx += edgeTable[i]["dNX"]
        jnx += edgeTable[j]["dNX"]
        iny += edgeTable[i]["dNY"]
        jny += edgeTable[j]["dNY"]
        inz += edgeTable[i]["dNZ"]
        jnz += edgeTable[j]["dNZ"]

        # if reached the bottom of an edge, switch out to next edge until reach bottom
        if (y >= edgeTable[i]["yEnd"] and y < lastFillLine):
            i = next
            ix = edgeTable[i]["xStart"]
            iz = edgeTable[i]["zStart"]
            inx = edgeTable[i]["nxStart"]
            iny = edgeTable[i]["nyStart"]
            inz = edgeTable[i]["nzStart"]
            next += 1
        if (y >= edgeTable[j]["yEnd"] and y < lastFillLine):
            j = next
            jx = edgeTable[j]["xStart"]
            jz = edgeTable[j]["zStart"]
            jnx = edgeTable[j]["nxStart"]
            jny = edgeTable[j]["nyStart"]
            jnz = edgeTable[j]["nzStart"]
            next += 1

# This function pre-computes the edge table to be used in Phong polygon shading (vertex normal values)
def computeEdgeTablePhongShade(polyDisplay,polyVertexNormals):
    # create a list of all the edges on the polygon, defined from point w/
    # smaller y value to point w/ larger y value, if the y values are equal
    # then don't add the edge - horizontal line
    edges = []
    for i in range(len(polyDisplay)):
        ptA = polyDisplay[i]
        ptB = polyDisplay[(i+1)%len(polyDisplay)]
        vA = polyVertexNormals[i]
        vB = polyVertexNormals[(i+1)%len(polyVertexNormals)]
        if (ptA[1] < ptB[1]):
            edges.append([ptA, ptB, vA, vB])
        elif (ptA[1] > ptB[1]):
            edges.append([ptB, ptA, vB, vA])
    
    # sort the edges in order of increasing yStart values
    # key is a function to sort by: a lambda function to find the yStart value of each edge
    edges.sort(key=(lambda edge : edge[0][1]))

    # create the final edge table, an array of dictionaries which store values and their var names
    # nx, ny, and nz are the surface normal values (x, y, and z)
    edgeTable = []
    for i in range(len(edges)):
        edgeTable.append({
            "edge": edges[i],
            "xStart": edges[i][0][0],
            "yStart": edges[i][0][1],
            "yEnd": edges[i][1][1],
            "dX": (edges[i][1][0]-edges[i][0][0])/(edges[i][1][1]-edges[i][0][1]),
            "zStart": edges[i][0][2],
            "dZ": (edges[i][1][2]-edges[i][0][2])/(edges[i][1][1]-edges[i][0][1]),
            "nxStart": edges[i][2][0],
            "dNX": (edges[i][3][0]-edges[i][2][0])/(edges[i][1][1]-edges[i][0][1]),
            "nyStart": edges[i][2][1],
            "dNY": (edges[i][3][1]-edges[i][2][1])/(edges[i][1][1]-edges[i][0][1]),
            "nzStart": edges[i][2][2],
            "dNZ": (edges[i][3][2]-edges[i][2][2])/(edges[i][1][1]-edges[i][0][1])
        })

    return edgeTable

# Generate a color hex code string from the illumination components
def triColorHexCode(ambient,diffuse,specular): 
    if diffuse < 0.00001:
        diffuse = 0
    if specular < 0.00001:
        specular = 0
    combinedColorCode = colorHexCode(ambient + diffuse + specular)
    specularColorCode = colorHexCode(specular)
    # keys 7-9 control the color of the cylinder (r,g,b) bools
    global rgb
    colorString = "#"
    for i in range(3):
        if rgb[i]:
            colorString += combinedColorCode
        else:
            colorString += specularColorCode
    return colorString

# Does the work of actually building the color string
def colorHexCode(intensity):
    hexString = str(hex(round(255*intensity)))
    if (hexString[0] == "-"): # illumination intensity should not be negative
        if (intensity < -0.1):
            print("Illumination intensity is negative. Setting to 00. Did you check for negative ndotl?")
        trimmedHexString = "00"
    else:
        trimmedHexString = hexString[2:] # get rid of "0x" at beginning of hex strings
        # convert single digit hex strings to two digit hex strings
        if len(trimmedHexString) == 1: trimmedHexString = "0" + trimmedHexString
        # we will use the green color component to display our monochrome illumination results
    return trimmedHexString

# Calculate a 3-D reflection vector, r, given surface normal, n, and lighting vetor, l
def reflect(n,l):
    r = []
    n = normalize(n)
    l = normalize(l)
    twoCosPhi = 2*(n[0]*l[0] + n[1]*l[1] + n[2]*l[2])
    if twoCosPhi > 0:
        for i in range(3):
            r.append(n[i] - (l[i]/twoCosPhi))
    elif twoCosPhi == 0:
        for i in range(3):
            r.append(-l[i])
    else: # twoCosPhi < 0
        for i in range(3):
            r.append(-n[i] + (l[i]/twoCosPhi))
    return normalize(r)

# **************************************************************************
# Everything below this point implements the interface
# For each function, delete all drawings, calculate the transformation to the 
# selected object, then draw the whole scene with the selected object red and bold.
def reset():
    w.delete(ALL)
    resetObject(currObject)
    drawScene(scene,currObject,renderMode)

def larger():
    w.delete(ALL)
    scale(scenePointClouds[currObject],1.1)
    drawScene(scene,currObject,renderMode)

def smaller():
    w.delete(ALL)
    scale(scenePointClouds[currObject],.9)
    drawScene(scene,currObject,renderMode)

def forward():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[0,0,5])
    drawScene(scene,currObject,renderMode)

def backward():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[0,0,-5])
    drawScene(scene,currObject,renderMode)

def left():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[-5,0,0])
    drawScene(scene,currObject,renderMode)

def right():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[5,0,0])
    drawScene(scene,currObject,renderMode)

def up():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[0,5,0])
    drawScene(scene,currObject,renderMode)

def down():
    w.delete(ALL)
    translate(scenePointClouds[currObject],[0,-5,0])
    drawScene(scene,currObject,renderMode)

def xPlus():
    w.delete(ALL)
    rotateX(scenePointClouds[currObject],5)
    drawScene(scene,currObject,renderMode)

def xMinus():
    w.delete(ALL)
    rotateX(scenePointClouds[currObject],-5)
    drawScene(scene,currObject,renderMode)

def yPlus():
    w.delete(ALL)
    rotateY(scenePointClouds[currObject],5)
    drawScene(scene,currObject,renderMode)

def yMinus():
    w.delete(ALL)
    rotateY(scenePointClouds[currObject],-5)
    drawScene(scene,currObject,renderMode)

def zPlus():
    w.delete(ALL)
    rotateZ(scenePointClouds[currObject],5)
    drawScene(scene,currObject,renderMode)

def zMinus():
    w.delete(ALL)
    rotateZ(scenePointClouds[currObject],-5)
    drawScene(scene,currObject,renderMode)

# This function reads which key the user pressed (left or right)
# and updates the currObject based on this input.
# Global had to be used because binding function in tkinter
# to read keyboard input cannot take any parameters.
def keyPress(key):
    global currObject, renderMode, rgb
    if (key.keycode == 39): # right key
        if (1 <= renderMode <= 3): # don't switch objects if in renderModes 4-6 (only cylinder)
            currObject = (currObject+1)%len(scene)
    elif (key.keycode == 37): # left key
        if (1 <= renderMode <= 3):
            currObject = (currObject-1)%len(scene)
    elif (key.keycode == 49): # 1 on keypad
        renderMode = 1
    elif (key.keycode == 50): # 2 on keypad
        renderMode = 2
    elif (key.keycode == 51): # 3 on keypad
        renderMode = 3
    elif (key.keycode == 52): # 4 on keypad
        renderMode = 4
        currObject = CYLINDER_SCENE_NO # only object in the scene, so make it selected
    elif (key.keycode == 53): # 5 on keypad
        renderMode = 5
        currObject = CYLINDER_SCENE_NO
    elif (key.keycode == 54): # 6 on keypad
        renderMode = 6
        currObject = CYLINDER_SCENE_NO
    elif (key.keycode == 55): # 7 on keypad - these 3 control the cylidner color in renderModes 4-6
        if (4 <= renderMode <= 6):
            rgb[0] = not(rgb[0])
    elif (key.keycode == 56): # 8 on keypad
        if (4 <= renderMode <= 6):
            rgb[1] = not(rgb[1])
    elif (key.keycode == 57): # 9 on keypad
        if (4 <= renderMode <= 6):
            rgb[2] = not(rgb[2])

    w.delete(ALL)
    drawScene(scene,currObject,renderMode)

root = Tk()
outerframe = Frame(root)
outerframe.pack()

w = Canvas(outerframe, width=CANVASWIDTH, height=CANVASHEIGHT)
w.delete(ALL)
drawScene(scene,currObject,renderMode)
w.pack()

controlpanel = Frame(outerframe)
controlpanel.pack()

resetcontrols = Frame(controlpanel, height=100, borderwidth=2, relief=RIDGE)
resetcontrols.pack(side=LEFT)

resetcontrolslabel = Label(resetcontrols, text="Reset")
resetcontrolslabel.pack()

resetButton = Button(resetcontrols, text="Reset", fg="green", command=reset)
resetButton.pack(side=LEFT)

scalecontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
scalecontrols.pack(side=LEFT)

scalecontrolslabel = Label(scalecontrols, text="Scale")
scalecontrolslabel.pack()

largerButton = Button(scalecontrols, text="Larger", command=larger)
largerButton.pack(side=LEFT)

smallerButton = Button(scalecontrols, text="Smaller", command=smaller)
smallerButton.pack(side=LEFT)

translatecontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
translatecontrols.pack(side=LEFT)

translatecontrolslabel = Label(translatecontrols, text="Translation")
translatecontrolslabel.pack()

forwardButton = Button(translatecontrols, text="FW", command=forward)
forwardButton.pack(side=LEFT)

backwardButton = Button(translatecontrols, text="BK", command=backward)
backwardButton.pack(side=LEFT)

leftButton = Button(translatecontrols, text="LF", command=left)
leftButton.pack(side=LEFT)

rightButton = Button(translatecontrols, text="RT", command=right)
rightButton.pack(side=LEFT)

upButton = Button(translatecontrols, text="UP", command=up)
upButton.pack(side=LEFT)

downButton = Button(translatecontrols, text="DN", command=down)
downButton.pack(side=LEFT)

rotationcontrols = Frame(controlpanel, borderwidth=2, relief=RIDGE)
rotationcontrols.pack(side=LEFT)

rotationcontrolslabel = Label(rotationcontrols, text="Rotation")
rotationcontrolslabel.pack()

xPlusButton = Button(rotationcontrols, text="X+", command=xPlus)
xPlusButton.pack(side=LEFT)

xMinusButton = Button(rotationcontrols, text="X-", command=xMinus)
xMinusButton.pack(side=LEFT)

yPlusButton = Button(rotationcontrols, text="Y+", command=yPlus)
yPlusButton.pack(side=LEFT)

yMinusButton = Button(rotationcontrols, text="Y-", command=yMinus)
yMinusButton.pack(side=LEFT)

zPlusButton = Button(rotationcontrols, text="Z+", command=zPlus)
zPlusButton.pack(side=LEFT)

zMinusButton = Button(rotationcontrols, text="Z-", command=zMinus)
zMinusButton.pack(side=LEFT)

# Reads for keyboard input to call the keyPress() function.
root.bind('<KeyPress>', keyPress)

root.mainloop()