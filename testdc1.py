#basic armature
from fury import window, actor
import numpy as np

scene = window.Scene()

# c = actor.line([np.random.rand(2, 3), np.random.rand(2, 3)])
p1 = [-1, 0, 0]
p2 = [1, 1, 0]
p3 = [-1, 1, 0]
p4 = [1, 0, 0]
p5 = [0.5, 0, 0]
p6 = [-0.5, 0, 0]
c = actor.line([[np.array(p1), np.array(p3)], [np.array(p4), np.array(p2)], 
                [np.array(p3), np.array(p2)], [np.array(p1), np.array(p6)], 
                [np.array(p5), np.array(p4)]], linewidth=10)

# dirs = np.array([[1, 1, 1]])
dirs = np.array([[0, 1, 0]])
colors = np.array([[1, 0, 0]])
centers = np.array([[0, 0, 0]])
scales = [[1]]

d = actor.sdf(centers=centers, directions=dirs, colors=colors,
              primitives=['torus'],
              scales=scales)

scene.add(c)
scene.add(d)
window.show(scene)