#only one ray

from fury import window, actor, utils, ui
import numpy as np
import itertools

time = 0
incre_time = 0.09
velocity = 1

mu1 = 2
mu2 = 1

scene = window.Scene()
scene.zoom(.4)
scene.set_camera(position=(10, 5, -20), focal_point=(1.0, 0.0, 0.0),
                 view_up=(0.0, 0.0, 0.0))
showm = window.ShowManager(scene,
                           size=(800, 600), reset_camera=True,
                           order_transparent=True)
showm.initialize()

yo1 = 10
zo1 = 10
xo1 = 10

dir_x = 1
dir_y = 1
dir_z = 1

dir1 = np.array([[dir_x, dir_y, dir_z]])
pts1 = np.array([[xo1, yo1, zo1]])

theta_i = np.arccos(np.dot(pts1[0], dir1[0]) / (np.sqrt(np.sum(pts1[0] ** 2)) * np.sqrt(np.sum(dir1[0] ** 2))))
theta_r = np.arcsin((mu2 / mu1) * np.sin(theta_i))
theta_c = np.arcsin(mu2 / mu1)

cos_theta_i = np.cos(theta_i)
cos_theta_r = np.cos(theta_r)

sin_theta_i = np.sin(theta_i)
sin_theta_r = np.sin(theta_r)

total_time1 = 2 * yo1 / (velocity * cos_theta_i)

color_particle = window.colors.red  


surface_points = np.array([[10, 0.02, 10],
                           [-10, 0.04, 10],
                           [-10, 0.03, -10],
                           [10, 0.01, -10]])


charge_actor1 = actor.point(np.array([[xo1, yo1, zo1]]), colors=(1, 0, 0))
surface = actor.surface(surface_points, colors=(0, 0, 0))

scene.add(charge_actor1)
scene.add(surface)

vertices1 = utils.vertices_from_actor(charge_actor1)
vcolors1 = utils.colors_from_actor(charge_actor1, 'colors')

no_vertices_per_point = len(vertices1)

initial_vertices1 = vertices1.copy() - np.repeat(pts1, no_vertices_per_point, axis=0)

counter = itertools.count()

end = 5000

coor_i1 = np.array([xo1, yo1, zo1])

def timer_callback1(_obj, _event):
    global pts, time, incre_time, coor_i1
    time += incre_time
    cnt = next(counter)

    if time < total_time1/2:
        y1 = yo1 - (velocity * time * cos_theta_i)
        z1 = zo1 - (velocity * time * cos_theta_i)
        x1 = xo1 - (velocity * time * cos_theta_i)

    elif (time >= total_time1/2):
        y1 = (velocity * time * cos_theta_i) - yo1
        z1 = zo1 - (velocity * time * cos_theta_i)
        x1 = xo1 - (velocity * time * cos_theta_i)

    pts1 = np.array([[x1, y1, z1]])

    vertices1[:] = initial_vertices1 + np.repeat(pts1, no_vertices_per_point, axis=0)

    if time < total_time1:
        utils.update_actor(charge_actor1)

        coor_f1 = np.array([x1, y1, z1])
        coors1 = np.array([coor_i1, coor_f1])
        coors1 = [coors1]
        line_actor = actor.line(coors1, window.colors.cyan, linewidth=3)
        scene.add(line_actor)
        coor_i1 = coor_f1

    showm.render()

    if cnt == end:
        showm.exit()


def timer_callback2(_obj, _event):
    global pts, time, incre_time, coor_i1
    time += incre_time
    cnt = next(counter)

    if time < total_time1/2:
        y1 = yo1 - (velocity * time * cos_theta_i)
        z1 = zo1 - (velocity * time * cos_theta_i)
        # x1 = xo1 - (velocity * time * cos_theta_i)
        x1 = -10

    elif (time >= total_time1/2):
        y1 = yo1 - (velocity * time * cos_theta_r)
        z1 = zo1 - (velocity * time * cos_theta_r)
        # x1 = xo1 - (velocity * time * cos_theta_r)
        x1 = -10
    pts1 = np.array([[x1, y1, z1]])

    vertices1[:] = initial_vertices1 + np.repeat(pts1, no_vertices_per_point, axis=0)

    if time < total_time1:
        utils.update_actor(charge_actor1)

        coor_f1 = np.array([x1, y1, z1])
        coors1 = np.array([coor_i1, coor_f1])
        coors1 = [coors1]
        line_actor = actor.line(coors1, window.colors.cyan, linewidth=3)
        scene.add(line_actor)
        coor_i1 = coor_f1

    showm.render()

    if cnt == end:
        showm.exit()


if theta_i > theta_c:
    showm.add_timer_callback(True, 15, timer_callback1)

else:
    showm.add_timer_callback(True, 15, timer_callback2)

interactive = True
if interactive:
    showm.start()