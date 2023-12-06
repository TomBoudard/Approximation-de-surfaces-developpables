import openmesh as om

def developabilityDetectFunction(mesh):
    ...
    return 0

#Thaïs
def computeMeshNormals(mesh):
    ...
    return 0

#Thaïs
def initHeap(mesh):
    """
    Take a mesh and returns a heap of the vertices according to their developability
    """
    return 0

def updateHeap(mesh, vertex):
    ...
    return 0

#Tom
def updateVertex(mesh, vertex):
    ...
    return 0

def updateNeighbour(mesh, vertex):
    ...
    return 0


def main():
    maxIter = 20
    nbIter = 0
    filename = "../Objects/sphere.off"
    mesh = om.read_trimesh(filename)

    developabilityDetectFunction(mesh)
    computeMeshNormals(mesh)
    heapVertex = initHeap(mesh)

    while(...):
        updateVertex(...)
        updateNeighbour(...)

        heapVertex = initHeap(mesh)

        nbIter += 1

    return 0

