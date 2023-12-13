import openmesh as om
import numpy as np
from tools import *
import heapq

def developabilityDetectFunction(mesh):
    developability = 0
    for vertex in mesh.vertices():
        if not mesh.is_boundary(vertex): #else developability = 0
            developability = 2 * np.pi
            neighbours = [vh for vh in mesh.vv(vertex)]
            for i in range(len(neighbours)):
                vh1 = neighbours[i]
                vh2 = neighbours[(i + 1) % len(neighbours)]
                vector1 = vector(mesh, vertex, vh1)
                vector2 = vector(mesh, vertex, vh2)
                angle = np.arccos(
                    np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
                developability -= angle
        developability = developability**2 #vraiment utile?
        mesh.set_vertex_property('developability', vertex, developability)
        #print(mesh.vertex_property('developability')[vertex.idx()])
    return 0


#Thaïs
def initDic(mesh):
    """
    Take a mesh and returns the vertex with the maximum developability
    """
    # We want a max heap instead of a min heap so we need to multiply the values by -1
    vertices = {vertex.idx() : mesh.vertex_property('developability')[vertex.idx()] for vertex in mesh.vertices()}
    vertex = max(vertices_heap, key = vertices_heap.get)
    #print(vertices_heap) debug
    return vertex

def updateDic(mesh, vertex, heap):
    ''' Update the developability of the vertex and its neighbours in the dictionary'''
    neighbours = [vh for vh in mesh.vv(vertex)]
    if not mesh.is_boundary(vertex): #else developability = 0
            developability = 2 * np.pi
            for i in range(len(neighbours)):
                vh1 = neighbours[i]
                vh2 = neighbours[(i + 1) % len(neighbours)]
                vector1 = vector(mesh, vertex, vh1)
                vector2 = vector(mesh, vertex, vh2)
                angle = np.arccos(
                    np.dot(vector1, vector2) / (np.sqrt(np.dot(vector1, vector1)) * np.sqrt(np.dot(vector2, vector2))))
                developability -= angle
        developability = developability**2 #vraiment utile?
        mesh.set_vertex_property('developability', vertex, developability)
    return 0

#Tom
def updateVertex(mesh, vertex):
    ...
    return 0

def updateNeighbour(mesh, vertex):
    ...
    return 0


def main():
    print("---------- Début main\n")
    maxIter = 20
    nbIter = 0
    filename = "../Objects/sphere.off"
    mesh = om.read_trimesh(filename)

    developabilityDetectFunction(mesh)
    mesh.update_vertex_normals()
    # print(mesh.vertex_normals()) #pour debug


    heapVertex = initHeap(mesh)

    # while(...):
    #     updateVertex(...)
    #     updateNeighbour(...)

    #     heapVertex = initHeap(mesh)

    #     nbIter += 1

    # return 0

if __name__ == "__main__":
    main() 