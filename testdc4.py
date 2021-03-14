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
p1 = [-1, 0, 0]
p2 = [1, 1, 0]
p3 = [-1, 1, 0]
p4 = [1, 0, 0]
p5 = [0.5, 0, 0]
p6 = [-0.5, 0, 0]
# charge_actor = actor.line([[np.array(p1), np.array(p3)], [np.array(p4), np.array(p2)], 
#                 [np.array(p3), np.array(p2)], [np.array(p1), np.array(p6)], 
#                 [np.array(p5), np.array(p4)]], linewidth=10)
# charge_actor1 = actor.line([np.random.rand(2, 3)], linewidth=1)
charge_actor1 = actor.line([[np.array(p1), np.array(p3)]], linewidth=1)


scene.add(charge_actor1)
# scene.add(charge_actor2)

vertices = utils.vertices_from_actor(charge_actor1)
vcolors = utils.colors_from_actor(charge_actor1, 'colors')
no_vertices_per_point = len(vertices)
print(no_vertices_per_point)
initial_vertices = vertices.copy() - np.repeat(pts, no_vertices_per_point, axis=0)

counter = itertools.count()

end = 5000


def timer_callback(_obj, _event):
    global pts, time, incre_time, coor_1
    time += incre_time
    cnt = next(counter)

    # x = initial_velocity*time + 0.5*acc*(time**2)
    # y = np.sin(10*angular_frq*time + phase_angle)
    # z = np.cos(10*angular_frq*time + phase_angle)
    
    # y = np.sin(10*angular_frq*time + phase_angle)
    z = np.cos(10*angular_frq*time + phase_angle)
    y = 0
    x = 0

    pts = np.array([[x, y, z]])

    vertices[:] = initial_vertices + np.repeat(pts, no_vertices_per_point, axis=0)

    utils.update_actor(charge_actor1)

    showm.render()

    if cnt == end:
        showm.exit()




# showm.add_timer_callback(True, 15, timer_callback)

interactive = True
if interactive:
    showm.start()