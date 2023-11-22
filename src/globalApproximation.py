import openmesh as om
import numpy as np

from tools import *


def developabilityDetectFunction(mesh):
    developability = 0
    for vertex in mesh.vertices():
        if not mesh.is_boundary(vertex):
            sum = 2*np.pi
            neighbours = [vh for vh in mesh.vv(vertex)]
            for i in range(len(neighbours)):
                vh1 = neighbours[i]
                vh2 = neighbours[(i+1)%len(neighbours)]
                vector1 = vector(mesh, vertex, vh1)
                vector2 = vector(mesh, vertex, vh2)
                angle = np.arccos(np.dot(vector1, vector2)/(np.sqrt(np.dot(vector1, vector1))*np.sqrt(np.dot(vector2, vector2))))
                sum -= angle
            developability += sum
        mesh.set_color(vertex, [developability, developability, developability, 1.])

    return developability, mesh


def main():
    filename = "../Objects/sphere.off "
    nbMaxIteration = 10

    mesh = om.read_trimesh(filename)

    developability, output = developabilityDetectFunction(mesh)
    nbIteration = 0

    print(developability)
    # while developability > 0 or nbIteration < nbMaxIteration:
    #     #TODO Change mesh
    #     developability = developabilityDetectFunction(mesh)
    #     nbIteration += 1


    om.write_mesh("Output.off", output, vertex_color=True)

    return 0


if __name__ == "__main__":
    main()