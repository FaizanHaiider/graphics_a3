#!/usr/bin/python

from objects import *
from math import *
import copy
from fractions import Fraction

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
    sphCoords = input_data[4].strip("\r\n").split(",")
    sphCoords[0], sphCoords[1], sphCoords[2] = int(sphCoords[0]), int(sphCoords[1]), int(sphCoords[2])
    return filename, imgSize, numShapes, shapeArr, triMesh, sphCoords
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
def createPolygonSurfaces(newPolygon, imgSize, origin):
    numSides = newPolygon.getNumSides()
    zVal = 200.0
    # create front face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, zVal)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # calculate z val for rear face
    zVal = -200.0

    # create rear face
    surfaceVertices = calcSurfaceVertices(imgSize, numSides, zVal)
    newPolygon.addSurface(surface(surfaceVertices, numSides))

    # create intermediate surfaces
    for i in range(numSides):
        newPolygon.addSurface(createInterSurfaces(newPolygon, i))
    return newPolygon
# end createPolygonSurfaces

# calculate surface vertices
def calcSurfaceVertices(imgSize, numSides, zVal):
    surfaceVertices = [0] * numSides
    for i in range(numSides):
        x = ((imgSize/2) * cos((2 * pi * i) / numSides) + (imgSize/2)) - (imgSize / 2)
        y = ((imgSize/2) * sin((2 * pi * i) / numSides) + (imgSize/2)) - (imgSize / 2)
        newVertex = vertex(x,y,zVal)
        surfaceVertices[i] = newVertex
    return surfaceVertices 
# end calcSurfaceVertices

# create surfaces between front/rear faces
def createInterSurfaces(newPolygon, i):
    numSides = newPolygon.getNumSides()
    tmpFront = newPolygon.getSurface(0)
    tmpRear = newPolygon.getSurface(1) 

    ii = i + 1
    if(ii == numSides):
        ii = 0
    vertices = [tmpRear.getVertex(i), tmpRear.getVertex(ii), tmpFront.getVertex(ii), tmpFront.getVertex(i)] 
    return surface(vertices, 4)
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
def rasterize(Matrix, poly, imgSize, viewpoint, origin):
    # modify all points so they are relative to the view point

    M, rM = createAugMatrix(viewpoint)
    det, rM = gauss(M, rM)
    newPoly = updatePolyPOV(rM, poly)

    for surface in newPoly.getSurfaceArr():
        verticesLen = surface.getVerticesLen()
        # regular edges
        for i in range(verticesLen):
            ii = i + 1
            if(ii == verticesLen):
                ii = 0
            v = surface.getVertex(i)
            vv = surface.getVertex(ii)
            BresenhamAlgo(v, vv, imgSize, Matrix, origin)
        # triangle mesh rasterize
        # for triangle in surface.getTriangles():
        #     v = triangle.getVertex(0)
        #     vv = triangle.getVertex(1)
        #     vvv = triangle.getVertex(2)
        #     BresenhamAlgo(v, vv, imgSize, Matrix)
        #     BresenhamAlgo(vv, vvv, imgSize, Matrix)
        #     BresenhamAlgo(vvv, v, imgSize, Matrix)
# end rasterize

# Bresenhams's Algo
def BresenhamAlgo(v, vv, imgSize, Matrix, origin):
    xs, ys = int(v.getx()) + imgSize/2, int(v.gety() + imgSize/2)
    xe, ye = int(vv.getx() + imgSize/2), int(vv.gety() + imgSize/2)

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

# find norm of vertex
def calcNorm(v):
    return sqrt(calcInnerProd(v,v))
# end calcNorm 

# scalar multiplication
def scalarMult(const, v):
    nv = vertex(v.getx() * const, v.gety() * const, v.getz() * const)
    return nv
# end scalarMult

# normalize unit vector
def normalize(v):
    # find norm
    norm = calcNorm(v)
    # constant
    const = (1/norm)
    # normalize
    return scalarMult(const, v)
# end normalize

# calc inner product
def calcInnerProd(v, vv):
    x = v.getx() * vv.getx()
    y = v.gety() * vv.gety()
    z = v.getz() * vv.getz()
    return x+y+z
# end calcInnerProd

# subtract two vectors
def subtrVector(v, vv):
    nvv = vertex(vv.getx() * (-1), vv.gety() * (-1), vv.getz() * (-1))
    fv = vertex(v.getx() + nvv.getx(), v.gety() + nvv.gety(), v.getz() + nvv.getz())
    return fv
# end subtrVector

# square function
def square(num):
    return num**2
# end square

# create viewPoint
def createViewPoint(sphCoords, origin, basis):
    # required variables
    r = sphCoords[0]
    elvation = radians(sphCoords[1])
    azi = radians(sphCoords[2])

    # calc C
    # calc x = r * cos(azi) * cos(elvation)
    x = r * sin(elvation) * cos(azi)
    # calc y = r * cos(azi) * sin(elvation)
    y = r * sin(elvation) * sin(azi)
    # calc z = r * cos(azi)
    z = r * cos(elvation)
    c = vertex(x,y,z)

    # calc N
    n = createVector(c, origin)
    n = normalize(n)

    # calc V
    k = basis.getk()
    kn = calcInnerProd(k, n)
    knn = scalarMult(kn, n)    
    vPrime = subtrVector(k, knn)
    v = normalize(vPrime)

    # calc U
    u = calcCrossProduct(n, v)

    viewpoint = viewPoint(c,u,v,n)
    return viewpoint
# end createViewPoint

# create aug matrix for gauss
def createAugMatrix(viewpoint):
    M = [[0 for x in range(3)] for y in range(3)]
    M[0][0] = viewpoint.getu().getx()
    M[1][0] = viewpoint.getu().gety()
    M[2][0] = viewpoint.getu().getz()

    M[0][1] = viewpoint.getv().getx()
    M[1][1] = viewpoint.getv().gety()
    M[2][1] = viewpoint.getv().getz()

    M[0][2] = viewpoint.getn().getx()
    M[1][2] = viewpoint.getn().gety()
    M[2][2] = viewpoint.getn().getz()

    resultM = [[0 for x in range(4)] for y in range(3)]

    resultM[0][3] = viewpoint.getc().getx() * -1
    resultM[1][3] = viewpoint.getc().gety() * -1
    resultM[2][3] = viewpoint.getc().getz() * -1 

    resultM[0][0], resultM[1][1], resultM[2][2] = 1, 1, 1
    resultM[0][1], resultM[0][2], resultM[1][0], resultM[1][2], resultM[2][0], resultM[2][1] = 0, 0, 0, 0, 0, 0
    
    return M, resultM
# end createAugMatrix

# update polygon coords with new POV coords
def updatePolyPOV(rM, poly):
    numSurfaces = poly.getNumSurfaces()
    surfaces = [0] * numSurfaces
    newPoly = polygon(numSurfaces, surfaces)

    for Surface in poly.getSurfaceArr():
        numVertices = Surface.getVerticesLen()
        vertices = [0] * numVertices
        index = 0
        for v in Surface.getVertices():
            nx, ny, nz = 0,0,0
            nx = (rM[0][0] * v.getx()) + (rM[0][1] * v.gety()) + (rM[0][2] * v.getz()) + rM[0][3]
            ny = (rM[1][0] * v.getx()) + (rM[1][1] * v.gety()) + (rM[1][2] * v.getz()) + rM[1][3]
            nz = (rM[2][0] * v.getx()) + (rM[2][1] * v.gety()) + (rM[2][2] * v.getz()) + rM[2][3]
            nv = vertex(nx, ny, nz)
            vertices[index] = nv
            index = index + 1
        newSurface = surface(vertices, numVertices)
        newPoly.addSurface(newSurface)
    return newPoly  
# end updatePolyPOV

# Gauss-Jordan Elimination
def gauss(a, b):
    a = copy.deepcopy(a)
    b = copy.deepcopy(b)
    n = len(a)
    p = len(b[0])
    det = 1
    for i in range(n - 1):
        k = i
        for j in range(i + 1, n):
            if abs(a[j][i]) > abs(a[k][i]):
                k = j
        if k != i:
            a[i], a[k] = a[k], a[i]
            b[i], b[k] = b[k], b[i]
            det = -det
 
        for j in range(i + 1, n):
            t = a[j][i]/a[i][i]
            for k in range(i + 1, n):
                a[j][k] -= t*a[i][k]
            for k in range(p):
                b[j][k] -= t*b[i][k]
 
    for i in range(n - 1, -1, -1):
        for j in range(i + 1, n):
            t = a[i][j]
            for k in range(p):
                b[i][k] -= t*b[j][k]
        t = 1/a[i][i]
        det *= a[i][i]
        for j in range(p):
            b[i][j] *= t
    return det, b
 
def zeromat(p, q):
    return [[0]*q for i in range(p)]
 
def matmul(a, b):
    n, p = len(a), len(a[0])
    p1, q = len(b), len(b[0])
    if p != p1:
        raise ValueError("Incompatible dimensions")
    c = zeromat(n, q)
    for i in range(n):
        for j in range(q):
                c[i][j] = sum(a[i][k]*b[k][j] for k in range(p))
    return c
 
def mapmat(f, a):
    return [list(map(f, v)) for v in a]
 
def ratmat(a):
    return mapmat(Fraction, a)
# end Gauss Jordon Elimination