import openmesh as om
import math 
import numpy as np
  
def dist(mesh,x1, x2):  
    x1 = mesh.point(x1)
    x2 = mesh.point(x2)
    xDiff = x1[0] - x2[0]  
    yDiff = x1[1] - x2[1]  
    zDiff = x1[2] - x2[2]  
    return np.sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff)

def compute_angles(mesh, A, B, C):  
    ''' Compute the angles of the given triangle'''
    # length of sides be a, b, c  
    a = dist(mesh, B,C)  
    b = dist(mesh, A,C)    
    c = dist(mesh, A,B)

    # Cosinus law  
    alpha = math.acos((b**2 + c**2 - a**2) / (2 * b * c))
    beta = math.acos((a**2 + c**2 - b**2) / (2 * a * c))
    gamma = math.acos((a**2 + b**2 - c**2) / (2 * a * b))

    return alpha,beta,gamma

def compute_area(mesh, face_vertices_dict, face_vertices_angles, face_id):
    ''' Compute the area of the face given'''
    A,B,C = face_vertices_dict[face_id]
    angle = face_vertices_angles[face_id][A.idx()]
    b = dist(mesh, A,C)
    c = dist(mesh, A,B)
    area = b * c * np.sin(angle) / 2
    return area
          

def add_angles(mesh):
    '''This function aim to store for every face of our mesh, the vertices associated and the angle associated to
    each verticies in the face'''
    face_vertices_dict = {} 
    face_vertices_angles = {}
    for face in mesh.faces(): #Iterate on the face of our mesh
        # print(f"Face {face.idx()}: ", end="")
        vertices = [v for v in mesh.fv(face)]  # Get vertex of the face
        face_vertices_dict[face.idx()] = vertices  # Store face ID as key and the verteices from the face as value

    for face in mesh.faces():
        vertices = face_vertices_dict[face.idx()] #get the 3 vertices our face
        a,b,c = compute_angles(mesh, vertices[0], vertices[1],vertices[2]) #comput the angles of our triangle
        vertices_angles = {vertices[0].idx(): a, vertices[1].idx() : b, vertices[2].idx(): c}
        face_vertices_angles[face.idx()] = vertices_angles #Store face ID as key and a list of the 3 [vertex ID, vertex's angle] as value    
    return face_vertices_dict, face_vertices_angles

def get_opposite_vertices(face_vertices_dict, face_id1, face_id2, vertex):
    '''this function return the vertex id of the vertex being in face1 and face2 but not being vertex_id 
    and also return the  two vertices ids of the two remaining vertices in the faces'''
    vertices_face1 = []
    vertices_face2 = []
    opp_vertices = []
    for i in range(3):
        vertex1 = face_vertices_dict[face_id1][i]
        vertex2 = face_vertices_dict[face_id2][i]
        if vertex1.idx() != vertex.idx():
            vertices_face1.append(vertex1.idx())
            if (vertex1.idx() in vertices_face2):
                opp_vertices.append(vertex1)
                vertices_face1.remove(vertex1.idx())
                vertices_face2.remove(vertex1.idx())
        if vertex2.idx() != vertex.idx(): 
            vertices_face2.append(vertex2.idx())
            if (vertex2.idx() in vertices_face1):
                opp_vertices.append(vertex2)
                vertices_face1.remove(vertex2.idx())
                vertices_face2.remove(vertex2.idx())
    opp_vertices.append(vertices_face1[0])
    opp_vertices.append(vertices_face2[0])
    return(opp_vertices)

def cotan(x):
    return np.cos(x)/np.sin(x)

def map_curvature_to_color(curvature_value):
    if curvature_value > 0:
        return (1.0, 0.0, 0.0, 1.0)  # Red color for positive curvature
    elif curvature_value < 0:
        return (0.0, 0.0, 1.0, 1.0)  # Blue color for negative curvature
    else:
        return (0.0, 1.0, 0.0, 1.0)  # Green color for zero curvature (or any default color)


def mean_curvature(mesh, face_vertices_dict, face_vertices_angles):
    #vertex_curvature = mesh.vertex_property("float")  
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
            #print(opp_vertices)
            #print(face_vertices_angles[face1])
            #print(face_vertices_angles[face1][opp_vertices[0]])
            cot_weight = cotan(face_vertices_angles[face1][opp_vertices[1]]) + cotan(face_vertices_angles[face2][opp_vertices[2]])
            pj = mesh.point(opp_vertices[0]) #get the coord. of the opposite vertex
            sum += cot_weight*(pj - mesh.point(vertex))
            face1_area = compute_area(mesh, face_vertices_dict, face_vertices_angles, face1)
            face2_area = compute_area(mesh, face_vertices_dict, face_vertices_angles, face2)
            Ai += face1_area + face2_area
        # print(face1_area)
        # print(face2_area)
        Ai *= 1/3
        # print("Ai = ", Ai)
        vertex_curvature = np.abs(np.linalg.norm(sum/(2*Ai))/2)
        #print("vertex_curvature =", vertex_curvature)
        # mesh.set_vertex_property('prop', vertex, vertex_curvature )
        color = map_curvature_to_color(vertex_curvature)
        mesh.set_color(vertex, (0.0, 1.0, 0.0, 1.0))
        # mesh.get_color(vertex)
    return mesh


######################################################################
def main2():   #Real main reading .off files
    print("---------- Début main\n")
    off_file = "../Objects/sphere.off"
    #à la main:
    #read_off_file(off_file)
    # En utilisnt les librairies:
    mesh = om.read_trimesh(off_file) #créer un mesh à partir du fichier .off
    #om.write_mesh('test.off', mesh) #Ecrit un fichier off à partir d'un mesh
    # TODO implémenter algorithme
    a, b = add_angles(mesh)
    mesh.request_face_colors()
    test = mean_curvature(mesh, a, b)
    om.write_mesh("test_sphere.ply", test)
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

    om.write_mesh("test_curvature.off", mesh) # write .off file
    #mean_curvature(mesh)
    print("---------- Fin main\n")

if __name__ == "__main__":
    main2()  # Call the main function when this file is executed directly