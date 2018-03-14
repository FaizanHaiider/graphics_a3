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
    
    def getVertices(self):
        return self.vertices
    def getVerticesLen(self):
        return self.numVertices
    def getVertex(self, index):
        return self.vertices[index]

# vertex/vector class
class vertex:
    def __init__(self, x, y, z):
        self.x = x
        self.y = y
        self.z = z
    
    def getx(self):
        return self.x
    def gety(self):
        return self.y
    def getz(self):
        return self.z
    def printVertex(self):
        print '({}, {}, {})'.format(self.x, self.y, self.z)
    