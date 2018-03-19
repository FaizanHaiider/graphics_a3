# polygon class
class polygon:
    def __init__(self, numSurfaces, surfaces):
        self.numSurfaces = numSurfaces
        self.surfaces = surfaces
        self.index = 0
        self.numSides = numSurfaces - 2
    
    def getNumSurfaces(self):
        return self.numSurfaces
    def getNumSides(self):
        return self.numSides
    def getSurfaceArr(self):
        return self.surfaces
    def getSurface(self, index):
        return self.surfaces[index]
    def addSurface(self, surface):
        self.surfaces[self.index] = surface
        self.index = self.index + 1

# surface class
class surface:
    def __init__(self, vertices, numVertices):
        self.vertices = vertices
        self.numVertices = numVertices
        self.triIndex = 0
        self.triangles = []
    
    def getVertices(self):
        return self.vertices
    def getVerticesLen(self):
        return self.numVertices
    def getVertex(self, index):
        return self.vertices[index]
    def getTriangles(self):
        return self.triangles
    def getNumTriangles(self):
        return self.triIndex
    def addTriangle(self, triangle):
        self.triangles.append(triangle)
        self.triIndex = self.triIndex + 1
    def printSurface(self):
        for vertex in self.vertices:
            vertex.printVertex()


# triangle class
class triangle:
    def __init__(self, vertices):
        self.vertices = vertices
        self.numVertices = 3
    
    def getVertices(self):
        return self.getVertices
    def getVertex(self, index):
        return self.vertices[index]
    def numVerticesLen(self):
        return self.numVertices

# vertex/vector class
class vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
        self.f = 1.0
    
    def getx(self):
        return self.x
    def gety(self):
        return self.y
    def getz(self):
        return self.z
    def setx(self, x):
        self.x = 0
        self.x = x
    def sety(self, y):
        self.y = 0
        self.y = y
    def setz(self, z):
        self.z = 0
        self.z = z
    def getf(self):
        return self.f
    def printVertex(self):
        print '({}, {}, {}, {})'.format(self.x, self.y, self.z, self.f)

# basis class
class Basis:
    def __init__(self, i, j, k):
        self.i = i
        self.j = j
        self.k = k
    
    def geti(self):
        return self.i
    def getj(self):
        return self.j
    def getk(self):
        return self.k

# view point class
class viewPoint:
    def __init__(self, c, u, v, n):
        self.c = c
        self.u = u
        self.v = v
        self.n = n
    
    def getc(self):
        return self.c
    def getu(self):
        return self.u
    def getv(self):
        return self.v
    def getn(self):
        return self.n
    def printViewPoint(self):
        self.c.printVertex()
        self.u.printVertex()
        self.v.printVertex()
        self.n.printVertex()