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

# update Matrix after rasterize()
def updateMatrix(Matrix, points, swapped):
    if(swapped):
        points.reverse()
    for point in points:
        Matrix[point[0]][point[1]] = 1
# end updateMatrix

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