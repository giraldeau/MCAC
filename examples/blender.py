from pathlib import Path

import bmesh
import bpy
import mathutils
import numpy as np

from pymcac import MCAC


def get_data():

    # The folder with all .h5 and .xmf files
    data_dir = Path("output_dir")

    # Read all data (~30s but the more data the longer)
    print("Reading")
    simu = MCAC(data_dir)
    simu.read_all_spheres()
    times = simu.times
    spheres = simu.spheres
    spheres = spheres.loc[times[-1]]
    spheres = spheres[spheres["Label"] == 0]

    norm = (
        max(
            spheres["Posx"].max() - spheres["Posx"].min(),
            spheres["Posy"].max() - spheres["Posy"].min(),
            spheres["Posz"].max() - spheres["Posz"].min(),
        )
        / 10
    )

    cx = spheres["Posx"].mean()
    cy = spheres["Posy"].mean()
    cz = spheres["Posz"].mean()

    positions = (
        np.stack((spheres["Posx"] - cx, spheres["Posy"] - cy, spheres["Posz"] - cz), axis=1) / norm
    )
    radius = np.array(spheres["Radius"]) / norm

    n = len(spheres["Radius"])
    colors = np.ones(n)

    data = [(colors[i], positions[i], radius[i]) for i in range(n)]
    return data


# create the materials used
# simple BI materials used
MatAu = bpy.data.materials.new("Mat.Au")
MatAu.diffuse_color = (0.8, 0.7, 0.2)
MatC = bpy.data.materials.new("Mat.C")
MatC.diffuse_color = (0.1, 0.1, 0.1)
MatH = bpy.data.materials.new("Mat.H")
MatH.diffuse_color = (0.8, 0.7, 0.2)

# create an empty mesh object and add it to the scene
sphereMesh = bpy.data.meshes.new("AllSpheres")
sphereObj = bpy.data.objects.new("AllSpheres", sphereMesh)
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
    mesh = bmesh.ops.create_uvsphere(
        bm, u_segments=8, v_segments=8, diameter=1.0, matrix=locMatrix * scaleMatrix
    )
    for v in mesh["verts"]:
        for f in v.link_faces:
            f.material_index = int(i[0])

bm.to_mesh(sphereMesh)
bm.free()
