import openmesh as om
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from tools import *
    
def mean_curvature(mesh, face_vertices_dict, face_vertices_angles):
    max_curvature = 0
    for vertex in mesh.vertices(): #iterate on the vertices of our mesh
        #print("New Vertex -----------------------------------------")
        sum = 0
        Ai = 0
        faces_ids = [face.idx() for face in mesh.vf(vertex)] #Store the ids of the faces adjacent to our vertex
        if len(faces_ids) == 1:
            #Vertex is only in one face -> Boundary vertex what shoud we do?
            #TODO
            continue
        for i in range(len(faces_ids)-1):
            #print(">>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>face n°", i)
            face1 = faces_ids[i]
            face2 = faces_ids[i+1]
            opp_vertices = get_opposite_vertices(face_vertices_dict, face1, face2, vertex)
            cot_weight = cotan(face_vertices_angles[face1][opp_vertices[1]]) + cotan(face_vertices_angles[face2][opp_vertices[2]])
            pj = mesh.point(opp_vertices[0]) #get the coord. of the opposite vertex
            sum += cot_weight*(pj - mesh.point(vertex))
            face1_area = compute_area(mesh, face_vertices_dict, face_vertices_angles, face1)
            face2_area = compute_area(mesh, face_vertices_dict, face_vertices_angles, face2)
            Ai += face1_area + face2_area
        Ai *= 1/3
        vertex_curvature = np.abs(np.linalg.norm(sum/(2*Ai))/2)
        if (vertex_curvature > max_curvature):
            max_curvature = vertex_curvature
        mesh.set_vertex_property('prop', vertex, vertex_curvature )
    for vertex in mesh.vertices():
        vertex_curvature = mesh.vertex_property('prop')[vertex.idx()]
        # print("vertex_curvature =", vertex_curvature)
        color = map_curvature_to_color(vertex_curvature, max_curvature)
        mesh.set_color(vertex, color) #à faire à posteriori après avoir trouvé max et min
    return mesh

######################################################################
def main2():   # More elaborate main reading and testing on .off files
    print("---------- Début main\n")
    off_file = "../../Objects/eight.off"
    # En utilisnt les librairies:
    mesh = om.read_trimesh(off_file) #créer un mesh à partir du fichier .off
    a, b = add_angles(mesh)
    test = mean_curvature(mesh, a, b)
    om.write_mesh("max_curvature.off", test, vertex_color = True)
    print("---------- Fin main\n")


def main():  # Dummy main with basic mesh to test the functions
    print("---------- Début main\n")
    mesh = om.TriMesh() #créer un mesh
    #initialize a mesh with 4 vertices
    vh0 = mesh.add_vertex([0, 1, 0])
    vh1 = mesh.add_vertex([1, 0, 0])
    vh2 = mesh.add_vertex([2, 1, 0])
    vh3 = mesh.add_vertex([0,-1, 0])
    # vh4 = mesh.add_vertex([2,-1, 0])

    # add a couple of faces to the mesh
    fh0 = mesh.add_face(vh0, vh1, vh2)
    #fh1 = mesh.add_face(vh1, vh3, vh4)
    fh2 = mesh.add_face(vh0, vh3, vh1)
    #vh_list = [vh2, vh1, vh4]
    #fh3 = mesh.add_face(vh_list)
    # TODO implémenter algorithme
    a, b = add_angles(mesh)
    mesh.request_face_colors() # pour ajouter des couleurs aux faces?
    new_mesh = mean_curvature(mesh, a, b)
    # mesh.get_color

    om.write_mesh("test_curvature.off", new_mesh, vertex_color = True) # write .off file
    #mean_curvature(mesh)
    print("---------- Fin main\n")

if __name__ == "__main__":
    main2()  # Call the main function when this file is executed directly