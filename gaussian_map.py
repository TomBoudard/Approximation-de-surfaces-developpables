import openmesh as om
import math 
import numpy as np
  
# returns square of distance b/w two openmesh vertices
def dist(mesh,x1, x2):  
    x1 = mesh.point(x1)
    x2 = mesh.point(x2)
    xDiff = x1[0] - x2[0]  
    yDiff = x1[1] - x2[1]  
    zDiff = x1[2] - x2[2]  
    return np.sqrt(xDiff * xDiff + yDiff * yDiff + zDiff * zDiff)

def comput_angles(mesh, A, B, C):  
      
    # length of sides be a, b, c  
    a = dist(mesh, B,C)  
    b = dist(mesh, A,C)    
    c = dist(mesh, A,B)
  
    # From Cosine law  
    alpha = math.acos((b**2 + c**2 - a**2) /
                         (2 * b * c));  
    betta = math.acos((a**2 + c**2 - b**2) / 
                         (2 * a * c));  
    gamma = math.acos((a**2 + b**2 - c**2) / 
                         (2 * a * b));  
  
    # Converting to degree  
    # alpha = alpha * 180 / math.pi;  
    # betta = betta * 180 / math.pi;  
    # gamma = gamma * 180 / math.pi;  
  
    return alpha,betta,gamma

def comput_area(mesh, face_vertices_dict, face_vertices_angles, face_id):
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
        vertices_ids = [v for v in mesh.fv(face)]  # Get vertex  of the face
        face_vertices_dict[face.idx()] = vertices_ids  # Store face ID as key and a list of the 3 [vertex ID, vertex's coord] as value

    for face in mesh.faces():
        vertices = face_vertices_dict[face.idx()] #get the 3 vertices id and coord of our face
        a,b,c = comput_angles(mesh, vertices[0], vertices[1],vertices[2]) #comput the angles of our triangle
        angles = [a, b, c] #store the 3 angles associated to our 3 vertices
        vertices_angles = {vertices[0].idx(): a, vertices[1].idx() : b, vertices[2].idx(): c}
        face_vertices_angles[face.idx()] = vertices_angles #Store face ID as key and a list of the 3 [vertex ID, vertex's angle] as value    
    return face_vertices_dict, face_vertices_angles

def get_opposite_vertices(face_vertices_dict, face_id1, face_id2, vertex):
    '''this function return the vertex id of the vertex also being in face1 and face2 but not being vertex_id 
    and also return the  two vertices ids of the vremaining vertices in the faces'''
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


def mean_curvature(mesh, face_vertices_dict, face_vertices_angles):
    #mesh is a openmesh mesh
    gaussian_curvature_values = []
    for vertex in mesh.vertices(): #iterate on the vertices of our mesh
        sum = 0
        Ai = 0
        faces_ids = [face.idx() for face in mesh.vf(vertex)] #Store the ids of the faces adjacent to our vertex
        for i in range(len(faces_ids)-1):
            #print("face n°", i)
            face1 = faces_ids[i]
            face2 = faces_ids[i+1]
            opp_vertices = get_opposite_vertices(face_vertices_dict, face1, face2, vertex)
            # print(opp_vertices)
           # print(face_vertices_angles[face1])
            #print(face_vertices_angles[face1][opp_vertices[0]])
            cot_weight = cotan(face_vertices_angles[face1][opp_vertices[1]]) + cotan(face_vertices_angles[face2][opp_vertices[2]])
            pj = mesh.point(opp_vertices[0]) 
            sum += cot_weight*(pj - mesh.point(vertex))
            face1_area = comput_area(mesh, face_vertices_dict, face_vertices_angles, face1)
            face2_area = comput_area(mesh, face_vertices_dict, face_vertices_angles, face2)
            print(face1_area)


############################################################################

# Assuming you have a triangle mesh defined by vertices and faces
# Replace these placeholders with your actual mesh data

vertices = [...]  # Your mesh vertices
faces = [...]     # Faces (triangles or polygons) defining the mesh

def compute_discrete_gaussian_curvature(mesh_vertices, mesh_faces):
    gaussian_curvature_values = []

    for vertex_idx, vertex in enumerate(mesh_vertices):
        adjacent_faces = np.where(np.any(mesh_faces == vertex_idx, axis=1))[0]

        vertex_gaussian_curvature = 0.0
        for face_idx in adjacent_faces:
            face = mesh_faces[face_idx]
            vertex_face_indices = np.where(face == vertex_idx)[0]

            if len(vertex_face_indices) == 0:
                continue

            opposite_vertex_indices = np.delete(face, vertex_face_indices)
            a = mesh_vertices[opposite_vertex_indices[0]] - vertex
            b = mesh_vertices[opposite_vertex_indices[1]] - vertex

            angle = np.arccos(np.dot(a, b) / (np.linalg.norm(a) * np.linalg.norm(b)))
            vertex_gaussian_curvature += angle

        # Calculate discrete Gaussian curvature
        vertex_gaussian_curvature = (2 * np.pi - vertex_gaussian_curvature) / np.linalg.norm(np.sum(mesh_vertices[mesh_faces[adjacent_faces]], axis=0) / len(adjacent_faces))

        gaussian_curvature_values.append(vertex_gaussian_curvature)

    return gaussian_curvature_values

# Calculate discrete Gaussian curvature values for your mesh
# discrete_gaussian_curvature_values = compute_discrete_gaussian_curvature(vertices, faces)
# print("Discrete Gaussian curvature values:", discrete_gaussian_curvature_values)


######################################################################
def main2():   #Vraie main qui utilisera des fichiers .off
    print("---------- Début main\n")
    off_file = "Objects/test.off"
    #à la main:
    #read_off_file(off_file)
    # En utilisnt les librairies:
    mesh = om.read_trimesh(off_file) #créer un mesh à partir du fichier .off
    #om.write_mesh('test.off', mesh) #Ecrit un fichier off à partir d'un mesh
    # TODO implémenter algorithme
    a, b = add_angles(mesh)
    mean_curvature(mesh, a, b)
    #mean_curvature(mesh)
    print("---------- Fin main\n")


def main():  # Dummy main avec mesh basique pour tester nos fonctions
    print("---------- Début main\n")
    mesh = om.TriMesh()

    # add a a couple of vertices to the mesh
    vh0 = mesh.add_vertex([0, 1, 0])
    vh1 = mesh.add_vertex([1, 0, 0])
    vh2 = mesh.add_vertex([2, 1, 0])
    vh3 = mesh.add_vertex([0,-1, 0])
    # vh4 = mesh.add_vertex([2,-1, 0])

    # add a couple of faces to the mesh
    fh0 = mesh.add_face(vh0, vh1, vh2)
    #fh1 = mesh.add_face(vh1, vh3, vh4)
    fh2 = mesh.add_face(vh0, vh3, vh1)
    # add another face to the mesh, this time using a list
    #vh_list = [vh2, vh1, vh4]
    #fh3 = mesh.add_face(vh_list)
    # TODO implémenter algorithme
    a, b = add_angles(mesh)
    mean_curvature(mesh, a, b)
    #mean_curvature(mesh)
    print("---------- Fin main\n")

if __name__ == "__main__":
    main()  # Call the main function when this file is executed directly