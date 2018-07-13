import bpy
import bmesh
import mathutils

import sys

sys.path.append(".")

from pathlib import Path
import numpy as np
import pandas as pd
from DLCA_Reader import DLCA

def get_data():

    # The folder with all .h5 and .xmf files
    data_dir = Path("output_dir")

    # Read all data (~30s but the more data the longer)
    print("Reading")
    simu = DLCA(data_dir)
    simu.read_all_spheres()
    times = simu.times
    Spheres = simu.Spheres
    Spheres = Spheres.loc[times[-1]]
    Spheres = Spheres[Spheres["Label"] == 0]

    norm = max(Spheres["Posx"].max() - Spheres["Posx"].min(),
               Spheres["Posy"].max() - Spheres["Posy"].min(),
               Spheres["Posz"].max() - Spheres["Posz"].min())/10

    Cx = Spheres["Posx"].mean()
    Cy = Spheres["Posy"].mean()
    Cz = Spheres["Posz"].mean()

    Positions = np.stack((Spheres["Posx"] - Cx,
                          Spheres["Posy"] - Cy,
                          Spheres["Posz"] - Cz), axis=1) / norm
    radius = np.array(Spheres["Radius"]) / norm

    n = len(Spheres["Radius"])
    colors = np.ones(n)

    data = [ (colors[i], Positions[i], radius[i]) for i in range(n)]
    return data






# create the materials used
# simple BI materials used
MatAu = bpy.data.materials.new('Mat.Au')
MatAu.diffuse_color = ((0.8,0.7,0.2))
MatC  = bpy.data.materials.new('Mat.C')
MatC.diffuse_color = ((0.1,0.1,0.1))
MatH  = bpy.data.materials.new('Mat.H')
MatH.diffuse_color = ((0.8,0.7,0.2))

# create an empty mesh object and add it to the scene
sphereMesh = bpy.data.meshes.new('AllSpheres')
sphereObj  = bpy.data.objects.new('AllSpheres', sphereMesh)
bpy.context.scene.objects.link(sphereObj)
bpy.context.scene.objects.active = sphereObj

# create slots for each material then add each material to the object
while len(sphereObj.material_slots) < 3:
    bpy.ops.object.material_slot_add()
sphereObj.material_slots[0].material = MatAu
sphereObj.material_slots[1].material = MatC
sphereObj.material_slots[2].material = MatH

# test data to be swapped for data in file
# type, location, scale
data = get_data()

bm = bmesh.new()
for i in data:
    locMatrix = mathutils.Matrix.Translation(i[1])
    scaleMatrix = mathutils.Matrix.Scale(i[2], 4)
    mesh = bmesh.ops.create_uvsphere(bm, u_segments=8, v_segments=8,
                    diameter=1.0, matrix=locMatrix*scaleMatrix)
    for v in mesh['verts']:
        for f in v.link_faces:
            f.material_index = int(i[0])

bm.to_mesh(sphereMesh)
bm.free()
