#!/usr/bin/python

from objects import *
from math import *

# read from input file
def readFileData(filename):
    input_data = []
    # read all lines
    with open(filename, "r") as f:
        input_data = f.readlines()
    
    # extract specific data
    filename = input_data[0].strip("\r\n")
    imgSize = int(input_data[1])
    numShapes, shapeArr = getShapeInfo(input_data[2].strip("\r\n"))
    triMesh = int(input_data[3].strip("\r\n"))
    return filename, imgSize, numShapes, shapeArr, triMesh
# end readFileData

# format shapes data from input file
def getShapeInfo(inputStr):
    # get number of shapes
    numShapes = inputStr.count(',') + 1
    shapeArr = [0] * numShapes

    if(numShapes) == 1:
        shapeArr[0] = int(inputStr)
        return numShapes, shapeArr
    
    inputStr = inputStr.split(",")
    for i in range(numShapes):
        shapeArr[i] = int(inputStr[i])
    return numShapes, shapeArr
# end getShapeInfo

# create surfaces for polygon
def createPolygonSurfaces(newPolygon, imgSize):
    numSides = newPolygon.getNumSides()
    # create front face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, 0)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # calculate z val for rear face
    tmpVecOne = createVector(surfaceVertices[0], surfaceVertices[1])
    tmpVecTwo = createVector(surfaceVertices[2], surfaceVertices[1])
    zVal = calcCrossProduct(tmpVecOne, tmpVecTwo).getz() / imgSize

    # create rear face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, zVal)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # create intermediate surfaces
    createInterSurfaces(newPolygon)
# end createPolygonSurfaces

# calculate surface vertices
def calcSurfaceVertices(imgSize, numSides, zVal):
    surfaceVertices = [0] * numSides
    for i in range(numSides):
        x = (imgSize/2) * cos((2 * pi * i) / numSides) + (imgSize/2)
        y = (imgSize/2) * sin((2 * pi * i) / numSides) + (imgSize/2)
        newVertex = vertex(x,y,zVal)
        surfaceVertices[i] = newVertex
    return surfaceVertices 
# end calcSurfaceVertices

# create surfaces between front/rear faces
def createInterSurfaces(newPolygon):
    numSides = newPolygon.getNumSides()
    tmpFront = newPolygon.getSurface(0)
    tmpRear = newPolygon.getSurface(1)

    for i in range(numSides):
        tmpI = i
        tmpII = i + 1
        if(tmpII == numSides):
            tmpII = 0
        vertices = [tmpRear.getVertex(tmpI), tmpRear.getVertex(tmpII), tmpFront.getVertex(tmpII), tmpFront.getVertex(tmpI)] 
        newPolygon.addSurface(surface(vertices, 4))
# end createInterSurfaces 

# create triangle mesh
def createTriangleMesh(polygon, imgSize, triMesh):
    numSurfaces = polygon.getNumSurfaces()
    numSides = polygon.getNumSides()
    startIndex = 0
    
    # simple tringle mesh for front/rear faces
    if(numSides > 4):
        createSimpleTriangleMesh(polygon.getSurface(0), imgSize)
        createSimpleTriangleMesh(polygon.getSurface(1), imgSize)
        startIndex = 2       
# end createTriangleMesh

# create simple triangle mesh
def createSimpleTriangleMesh(surface, imgSize):
    # calc center
    z = surface.getVertex(0).getz()
    center = vertex(int(imgSize/2), int(imgSize/2), z)

    # simple triangle mesh
    for i in range(surface.getVerticesLen()-1):
        vertices = [0] * 3
        v = surface.getVertex(i)
        vv = surface.getVertex(i+1)
        vertices[0], vertices[1], vertices[2] = v, vv, center
        ntriangle = triangle(vertices)
        surface.addTriangle(ntriangle)
# end createSimpleTriangleMesh

# rasterize polygon
def rasterize(Matrix, poly, imgSize):
    #
    for surface in poly.getSurfaceArr():
        verticesLen = surface.getVerticesLen()
        # regular edges
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            v = surface.getVertex(i)
            vv = surface.getVertex(ii)
            BresenhamAlgo(v, vv, imgSize, Matrix)
        # triangle mesh rasterize
        for triangle in surface.getTriangles():
            v = triangle.getVertex(0)
            vv = triangle.getVertex(1)
            vvv = triangle.getVertex(2)
            BresenhamAlgo(v, vv, imgSize, Matrix)
            BresenhamAlgo(vv, vvv, imgSize, Matrix)
            BresenhamAlgo(vvv, v, imgSize, Matrix)
# end rasterize

# Bresenhams's Algo
def BresenhamAlgo(v, vv, imgSize, Matrix):
    xs, ys = int(v.getx()), int(v.gety())
    xe, ye = int(vv.getx()), int(vv.gety())

    dx = xe - xs
    dy = ye - ys
    isSteep = abs(dy) > abs(dx)

    if(isSteep):
        xs, ys = ys, xs
        xe, ye = ye, xe
    swapped = 0
    if(xs > xe):
        xs, xe = xe, xs
        ys, ye = ye, ys
        swapped = 1
    
    dx = xe - xs
    dy = ye - ys
    error = int(dx/2.0)
    yincr = 1 if(ys < ye) else -1

    y = ys
    points = []
    for x in range(xs, xe + 1):
        tx = (y if isSteep else x)
        ty = (x if isSteep else y)
        tx = (tx if tx < imgSize else tx-1)
        ty = (ty if ty < imgSize else ty-1)
        coord = (tx, ty)
        points.append(coord)
        error  = error - abs(dy)
        if(error < 0):
            y = y + yincr
            error = error + dx
    updateMatrix(Matrix, points, swapped)    
# end BresenhamAlgo

# update Matrix after rasterize()
def updateMatrix(Matrix, points, swapped):
    if(swapped):
        points.reverse()
    for point in points:
        Matrix[point[0]][point[1]] = 1
# end updateMatrix

# write Matrix to pbm file
def writePPM(Matrix, filename, imgSize):
    # create output file
    with open(filename, "w+") as fp:
        for i in Matrix:
            fp.write(' '.join(map(str, i)))

    # read, format, write file data
    with open(filename, "r+") as fp:
        lines = fp.read()
    lines = lines.replace("00", "0 0")
    lines = lines.replace("01", "0 1")
    lines = lines.replace("10", "1 0")
    lines = lines.replace("11", "1 1")
    with open(filename, "w+") as fp:
        fp.write("P1\n")
        fp.write('{} {}\n'.format(imgSize, imgSize))
        fp.write(lines)
# end writePPM

# calculate cross product two vectors
def calcCrossProduct(vOne, vTwo):
	x = (vOne.gety() * vTwo.getz()) - (vTwo.gety() * vOne.getz())
	y = (vOne.getz() * vTwo.getx()) - (vTwo.getz() * vOne.getx())
	z = (vOne.getx() * vTwo.gety()) - (vTwo.getx() * vOne.gety())
	newVertex = vertex(x,y,z)
	return newVertex
# end calcCrossProd

# create vector using two vertices
def createVector(vOne, vTwo):
    x = vTwo.getx() - vOne.getx()
    y = vTwo.gety() - vOne.gety()
    z = vTwo.getz() - vOne.getz()

    return vertex(x,y,z)
# end createVector

# square function
def square(num):
    return num**2