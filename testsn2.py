##rotaion
from fury import window, actor, utils, ui
import numpy as np
import itertools

time = 0
incre_time = 0.09
velocity = 1

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



yo1 = 10
zo1 = 10
xo1 = 10

total_time1 = 2 * yo1 / (velocity / np.sqrt(3))

yo2 = 10
zo2 = 1
xo2 = 1

color_particle = window.colors.red  # color of particle can be manipulated
pts1 = np.array([[xo1, yo1, zo1]])
pts2 = np.array([[xo2, yo2, zo2]])
# pts2 = np.array([[x2, y2, z2]])

# line1 = [np.array([[-2.0, 0.0, 0.0], [2.0, 0.0, 0.0]])]
# line2 = [np.array([[-2.0, 1.0, 0.0], [2.0, 1.0, 0.0]])]

charge_actor1 = actor.point(np.array([[xo1, yo1, zo1]]), colors=(1, 0, 0))
charge_actor2 = actor.point(np.array([[xo2, yo2, zo2]]), colors=(1, 0, 0))

scene.add(charge_actor1)
scene.add(charge_actor2)

vertices1 = utils.vertices_from_actor(charge_actor1)
vcolors1 = utils.colors_from_actor(charge_actor1, 'colors')

vertices2 = utils.vertices_from_actor(charge_actor2)
vcolors2 = utils.colors_from_actor(charge_actor2, 'colors')

no_vertices_per_point = len(vertices1)

initial_vertices1 = vertices1.copy() - np.repeat(pts1, no_vertices_per_point, axis=0)
initial_vertices2 = vertices2.copy() - np.repeat(pts2, no_vertices_per_point, axis=0)

counter = itertools.count()

end = 5000

coor_i1 = np.array([xo1, yo1, zo1])
coor_i2 = np.array([xo2, yo2, zo2])

def timer_callback(_obj, _event):
    global pts, time, incre_time, coor_i1, coor_i2
    time += incre_time
    cnt = next(counter)

    # x = initial_velocity*time + 0.5*acc*(time**2)
    # y = np.sin(10*angular_frq*time + phase_angle)
    # z = np.cos(10*angular_frq*time + phase_angle)
    
    if time < total_time1/2:
        y1 = yo1 - (velocity * time)/(np.sqrt(3))
        z1 = zo1 - (velocity * time)/(np.sqrt(3))
        x1 = xo1 - (velocity * time)/(np.sqrt(3))
    
    elif (time >= total_time1/2):
        y1 = (velocity * time)/(np.sqrt(3)) - yo1
        z1 = zo1 - (velocity * time)/(np.sqrt(3))
        x1 = xo1 - (velocity * time)/(np.sqrt(3))

    if time < total_time1/2:
        y2 = yo2 - (velocity * time)/(np.sqrt(3))
        z2 = zo2 - (velocity * time)/(np.sqrt(3))
        x2 = xo2 - (velocity * time)/(np.sqrt(3))
    
    elif (time >= total_time1/2):
        y2 = (velocity * time)/(np.sqrt(3)) - yo2
        z2 = zo2 - (velocity * time)/(np.sqrt(3))
        x2 = xo2 - (velocity * time)/(np.sqrt(3))

    # y2 = np.sin(10*angular_frq2*time + phase_angle2)
    # z2 = np.cos(10*angular_frq2*time + phase_angle2)
    # x2 = 0

    pts1 = np.array([[x1, y1, z1]])
    pts2 = np.array([[x2, y2, z2]])

    vertices1[:] = initial_vertices1 + np.repeat(pts1, no_vertices_per_point, axis=0)
    vertices2[:] = initial_vertices2 + np.repeat(pts2, no_vertices_per_point, axis=0)
    
    if time < total_time1:
        utils.update_actor(charge_actor1)
        utils.update_actor(charge_actor2)

    # utils.update_actor(charge_actor2)
        coor_f1 = np.array([x1, y1, z1])
        coors1 = np.array([coor_i1, coor_f1])
        coors1 = [coors1]
        line_actor = actor.line(coors1, window.colors.cyan, linewidth=3)
        scene.add(line_actor)
        coor_i1 = coor_f1

        coor_f2 = np.array([x2, y2, z2])
        coors2 = np.array([coor_i2, coor_f2])
        coors2 = [coors2]
        line_actor = actor.line(coors2, window.colors.cyan, linewidth=3)
        scene.add(line_actor)
        coor_i2 = coor_f2
    
    
    showm.render()

    if cnt == end:
        showm.exit()




showm.add_timer_callback(True, 15, timer_callback)

interactive = True
if interactive:
    showm.start()