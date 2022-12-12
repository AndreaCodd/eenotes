# Import modules:
import gmsh
import sys
 
# Initialize gmsh:
gmsh.initialize()

# Points
lc=0.1;
point1 = gmsh.model.geo.add_point(0, 0, 0, lc)
point2 = gmsh.model.geo.add_point(1, 0, 0, lc)
point3 = gmsh.model.geo.add_point(1, 1, 0, 0.2*lc)
point4 = gmsh.model.geo.add_point(0, 1, 0, 0.1*lc)

# Lines
line1 = gmsh.model.geo.add_line(point1, point2)
line2 = gmsh.model.geo.add_line(point2, point3)
line3 = gmsh.model.geo.add_line(point3, point4)
line4 = gmsh.model.geo.add_line(point4, point1)

# Surface
face1 = gmsh.model.geo.add_curve_loop([line1, line2, line3, line4])
gmsh.model.geo.add_plane_surface([face1])

# Make model, mesh, and write msh file
gmsh.model.geo.synchronize()
gmsh.model.mesh.generate()
gmsh.write("pythonGmshSquare.msh")
 

