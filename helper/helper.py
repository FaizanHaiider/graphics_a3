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
    triMesh = input_data[3].strip("\r\n")
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
    # createInterSurfaces(newPolygon)
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

# Bresenhams's Algo
def rasterize(Matrix, poly, imgSize):
    for surface in ploy.getSurfaceArr():
        verticesLen = surface.getVerticesLen()
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            v = surface.getVertex(i)
            vv = surface.getVertex(ii)

            xs, ys = v.getx(), v.gety()
            xe, ye = vv.getx(), vv.gety()

            dx = xe - xs
            dy = ye - ys
            isSteep = abs(dy) - abs(dx)

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
                coord = (y, x) if isSteep else (x, y)
                points.append(coord)
                error  = error - abs(dy)
                if(error < 0):
                    y = y + ystep
                    error = error + dx
            updateMatrix(points, swapped)    
# end rasterize

# update Matrix after rasterize()
def updateMatrix(points, swapped):
    if(swapped):
        points.reverse()
    for point in points:
        Matrix[point[0], point[1]] = 1
# end updateMatrix

# write Matrix to ppm file
def writePPM(Matrix, filename, imgSize):
    # create output file
    with open(filename, "w+") as fp:
        fp.write("P1")
        fp.write('{} {}'.format(imgSize, imgSize))
        for i in Matrix:
            fp.write(' '.join(map(str, i)))
# end writePPM

# calculate cross product two vectors
def calcCrossProduct(vOne, vTwo):
	x = (vOne.gety() * vTwo.getz()) - (vTwo.gety() * vOne.getz())
	y = (vOne.getz() * vTwo.getx()) - (vTwo.getz() * vOne.getx())
	z = (vOne.getx() * vTwo.gety()) - (vTwo.getx() * vOne.gety())
	newVertex = vertex(x,y,z)
	return newVertex
# end calcCrossProd

# create vector
def createVector(vOne, vTwo):
    x = vTwo.getx() - vOne.getx()
    y = vTwo.gety() - vOne.gety()
    z = vTwo.getz() - vOne.getz()

    return vertex(x,y,z)
# end createVector

# square function
def square(num):
    return num**2
