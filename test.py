import numpy as np
from fury import window, actor


dirs = np.random.rand(1, 3)
colors = np.random.rand(1, 3) * 255
centers = np.array([[0, 0, 0]])
scales = np.random.rand(1, 1)


# sdfactor = actor.sdf(centers=centers, directions=dirs, colors=colors,
#                      primitives=['sphere', 'torus', 'ellipsoid'],
#                      scales=scales)

sdfactor = actor.sdf(centers=centers, directions=dirs, colors=colors,
                     primitives=['se'],
                     scales=scales)


scene = window.Scene()
scene.background((1.0, 0.8, 0.8))
scene.add(sdfactor)


current_size = (1024, 720)
showm = window.ShowManager(scene, size=current_size,
                           title="Visualize SDF Actor")

interactive = True

if interactive:
    showm.start()

window.record(scene, out_path='viz_sdfactor.png', size=current_size)