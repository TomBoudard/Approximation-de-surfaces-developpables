import openmesh as om
import numpy as np
# import matplotlib.pyplot as plt

# Load a mesh (replace 'path/to/your/mesh.obj' with your mesh file)
off_file = "Objects/test.off"
mesh = om.read_trimesh(off_file)

# Access vertices and faces
vertex_array = np.array([[mesh.point(v)[0], mesh.point(v)[1], mesh.point(v)[2]] for v in mesh.vertices()])
face_array = np.array([[v.idx() for v in mesh.fv(fh)] for fh in mesh.faces()])

# Display the mesh using Matplotlib
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')

# # Plotting vertices
# ax.scatter(vertex_array[:, 0], vertex_array[:, 1], vertex_array[:, 2], c='blue', marker='.')

# # Plotting faces
# for face in face_array:
#     face_vertices = vertex_array[face]
#     x, y, z = face_vertices[:, 0], face_vertices[:, 1], face_vertices[:, 2]
#     ax.plot(x, y, z, color='black')

# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')

# plt.show()
