##rotaion
from fury import window, actor, utils, ui
import numpy as np
import itertools

radius_particle = 0.08
initial_velocity = 0.09
acc = 0.004
time = 0
incre_time = 0.09
# angular_frq = 0.1
angular_frq = 0.1
phase_angle = 0.002

scene = window.Scene()
scene.zoom(1.2)
scene.set_camera(position=(10, 12.5, 19), focal_point=(3.0, 0.0, 0.0),
                 view_up=(0.0, 0.0, 0.0))
showm = window.ShowManager(scene,
                           size=(800, 600), reset_camera=True,
                           order_transparent=True)
showm.initialize()


x = initial_velocity*time + 0.5*acc*(time**2)
y = np.sin(angular_frq*time + phase_angle)
z = np.cos(angular_frq*time + phase_angle)


color_particle = window.colors.red  # color of particle can be manipulated
pts = np.array([[x, y, z]])
# charge_actor = actor.point(pts, color_particle, point_radius=radius_particle)

# ini = [np.array([[1.0, 0.0, 0.0], [0.0, 0.0, 0.0]])]
ini1 = [np.array([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])]
ini2 = [np.array([[-2.0, 1.0, 0.0], [2.0, 1.0, 0.0]])]
# ini2 = [np.array([[0.0, 1.0, 0.0], [0.0, 0.0, 0.0]])]
# box_center1 = np.array([[0.0, 0.0, 5.0]])
box_center1 = np.array([[0.0, 0.0, 2.0]])
box_dirs1 = np.array([[1.0, 1.0, 0.0]])

box_center2 = np.array([[0.0, 0.0, -5.0]])
box_dirs2 = np.array([[1.0, 1.0, 0.0]])

charge_actor1 = actor.line(ini1, linewidth=1)
charge_actor2 = actor.line(ini2, linewidth=1)
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
initial_vertices1 = vertices1.copy() - np.repeat(pts, no_vertices_per_point, axis=0)
# initial_vertices2 = vertices2.copy() - np.repeat(pts, no_vertices_per_point, axis=0)
initial_vertices2 = vertices2.copy() - np.repeat(pts, 2, axis=0)

counter = itertools.count()

end = 5000


def timer_callback(_obj, _event):
    global pts, time, incre_time, coor_1
    time += incre_time
    cnt = next(counter)

    # x = initial_velocity*time + 0.5*acc*(time**2)
    # y = np.sin(10*angular_frq*time + phase_angle)
    # z = np.cos(10*angular_frq*time + phase_angle)
    
    y = np.sin(10*angular_frq*time + phase_angle)
    z = np.cos(10*angular_frq*time + phase_angle)
    # y = 0
    x = 0

    pts = np.array([[x, y, z]])

    vertices1[:] = initial_vertices1 + np.repeat(pts, no_vertices_per_point, axis=0)
    vertices2[:] = initial_vertices2 + np.repeat(pts, no_vertices_per_point, axis=0)

    utils.update_actor(charge_actor1)
    # utils.update_actor(charge_actor2)

    showm.render()

    if cnt == end:
        showm.exit()


showm.add_timer_callback(True, 15, timer_callback)

interactive = True
if interactive:
    showm.start()