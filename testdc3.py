##rotaion
from fury import window, actor, utils, ui
import numpy as np
import itertools

time = 0
incre_time = 0.09

angular_frq1 = 0.01
angular_frq2 = 0.01

phase_angle1 = 0.0
phase_angle2 = 3.14

scene = window.Scene()
scene.zoom(1.2)
scene.set_camera(position=(10, 12.5, 19), focal_point=(3.0, 0.0, 0.0),
                 view_up=(0.0, 0.0, 0.0))
showm = window.ShowManager(scene,
                           size=(800, 600), reset_camera=True,
                           order_transparent=True)
showm.initialize()

# y1 = np.sin(10*angular_frq1*time + phase_angle1)
# z1 = np.cos(10*angular_frq1*time + phase_angle1)
# x1 = 0

# y2 = np.sin(10*angular_frq2*time + phase_angle2)
# z2 = np.cos(10*angular_frq2*time + phase_angle2)
# x2 = 0


y1 = 0
z1 = 0
x1 = 0

y2 = 0
z2 = 0
x2 = 0

color_particle = window.colors.red  # color of particle can be manipulated
pts1 = np.array([[x1, y1, z1]])
pts2 = np.array([[x2, y2, z2]])

# line1 = [np.array([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])]
# line2 = [np.array([[-2.0, 1.0, 0.0], [2.0, 1.0, 0.0]])]


line1 = [np.array([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])]
line2 = [np.array([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])]


box_center1 = np.array([[0.0, 0.0, 3.0]])
box_dirs1 = np.array([[1.0, 1.0, 0.0]])

box_center2 = np.array([[0.0, 0.0, -3.0]])
box_dirs2 = np.array([[1.0, 1.0, 0.0]])

charge_actor1 = actor.line(line1, linewidth=1)
charge_actor2 = actor.line(line2, linewidth=1)
box1 = actor.box(box_center1, box_dirs1, (1, 0, 0))
box2 = actor.box(box_center2, box_dirs1, (1, 0, 0))


scene.add(charge_actor1)
scene.add(charge_actor2)
scene.add(box1)
scene.add(box2)

vertices1 = utils.vertices_from_actor(charge_actor1)
vcolors1 = utils.colors_from_actor(charge_actor1, 'colors')

vertices2 = utils.vertices_from_actor(charge_actor2)
vcolors2 = utils.colors_from_actor(charge_actor2, 'colors')

no_vertices_per_point = len(vertices1)
# print(no_vertices_per_point)
initial_vertices1 = vertices1.copy() - np.repeat(pts1, no_vertices_per_point, axis=0)
initial_vertices2 = vertices2.copy() - np.repeat(pts2, no_vertices_per_point, axis=0)

counter = itertools.count()

end = 5000


def timer_callback(_obj, _event):
    global pts, time, incre_time, coor_1
    time += incre_time
    cnt = next(counter)

    # x = initial_velocity*time + 0.5*acc*(time**2)
    # y = np.sin(10*angular_frq*time + phase_angle)
    # z = np.cos(10*angular_frq*time + phase_angle)
    
    y1 = np.sin(10*angular_frq1*time + phase_angle1)
    z1 = np.cos(10*angular_frq1*time + phase_angle1)
    x1 = 0

    y2 = np.sin(10*angular_frq2*time + phase_angle2)
    z2 = np.cos(10*angular_frq2*time + phase_angle2)
    x2 = 0

    pts1 = np.array([[x1, y1, z1]])
    pts2 = np.array([[x2, y2, z2]])

    vertices1[:] = initial_vertices1 + np.repeat(pts1, no_vertices_per_point, axis=0)
    vertices2[:] = initial_vertices2 + np.repeat(pts2, no_vertices_per_point, axis=0)

    utils.update_actor(charge_actor1)
    utils.update_actor(charge_actor2)

    showm.render()

    if cnt == end:
        showm.exit()




showm.add_timer_callback(True, 15, timer_callback)

interactive = True
if interactive:
    showm.start()