#/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np
import math

# inputed param
total_step   = 100000
save_step    = 1000
core_patch_set_num = 100
check_times  = 100
v0           = 1.0
tolerance    = 1e-6
# box size
x_edge       = 10.0
y_edge       = 2.0
z_edge       = 2.0

saved_step_num = total_step / save_step

traj         = open("data/test_g-BAOABIntegrator_position.xyz")
cores_traj   = list(filter(lambda x: x[0:4] == 'CORE',  traj))
traj.seek(0)
patches_traj = list(filter(lambda x: x[0:5] == 'PATCH', traj))

dist_log = []

for check_itr in range(1, check_times):
    check_step     = np.random.randint(saved_step_num)
    check_particle = np.random.randint(core_patch_set_num)
    checked_idx    = check_step * core_patch_set_num + check_particle
    core_vec    = np.array([ float(cores_traj  [checked_idx].split()[i]) for i in range(1,4) ])
    patch_vec   = np.array([ float(patches_traj[checked_idx].split()[i]) for i in range(1,4) ])
    dist_vec    = core_vec - patch_vec
    # fix periodic boundary condition effect
    if x_edge / 2.0 < math.fabs(dist_vec[0]):
        dist_vec[0] -= x_edge * np.sign(dist_vec[0])
    if y_edge / 2.0 < math.fabs(dist_vec[1]):
        dist_vec[1] -= y_edge * np.sign(dist_vec[1])
    if z_edge / 2.0 < math.fabs(dist_vec[2]):
        dist_vec[2] -= z_edge * np.sign(dist_vec[2])
    distance = np.linalg.norm(dist_vec)
    if(distance < v0 - tolerance or v0 + tolerance < v0):
        print(
              """
              constraint condition was broken at particle {particle} in step {step}
              distance is {distance:.4f}
              core coordinate is {core_x:.4f} {core_y:.4f} {core_z:.4f}
              patch coordinate is {patch_x:.4f} {patch_y:.4f} {patch_z:.4f}
              """.format(particle=check_particle, step=check_step, distance=distance,
                         core_x=core_vec[0],   core_y=core_vec[1],   core_z=core_vec[2],
                         patch_x=patch_vec[0], patch_y=patch_vec[1], patch_z=patch_vec[2])
              )
