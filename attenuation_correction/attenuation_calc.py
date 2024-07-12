import numpy as np
import matplotlib.pyplot as plt


def sq_to_circ(l):

    area = l**2
    r = np.sqrt(area/np.pi)
    return r


def point_calc(det_dist, det_R, N):

    counts = 0
    paths = np.zeros(shape=(N,3))
    tally = [0,0,0,0,0,0]

    for i in range(N):

        # Generate start position
        # start = np.array([side_l*np.random.random() - side_l/2, side_l*np.random.random() - side_l/2, thickness*np.random.random()])
        start = np.array([0, 0, 0])
        # Generate a random directional vector
        theta = np.random.uniform(0, 2 * np.pi)  # azimuthal angle
        # phi = np.arccos(1 - 2 * np.random.uniform(0, 1))  # polar angle
        phi = np.arccos(-np.random.uniform(0, 1))  # polar angle for downward hemisphere
        vec = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])

        # A photon travelling upwards will not be able to enter the detector
        if vec[2] > 0:
            print('fail')
            pass
        # Otherwise, check whether it will be possible to enter the detector
        else:
            det_vec = start + np.abs((det_dist+start[2])/vec[2])*vec # vector to detector plane
            if ((det_vec[0]**2 + det_vec[1]**2) > det_R**2):
                # print('miss', det_vec)
                pass
            else:
                # print('hit', det_vec)
                bottom_pos = start + np.abs(start[2]/vec[2])*vec
                path = bottom_pos - start
                distance = np.sqrt(path[0]**2 + path[1]**2 + path[2]**2)
                counts += 1

    return counts


def foil_calc_sq(side_l, thickness, det_dist, det_R, N, attenuate = False, mat='dy', E=100.10595):

    counts = 0
    paths = np.zeros(shape=(N,3))
    tally = [0,0,0,0,0,0]
    xcom=np.loadtxt(f'xcom_{mat}.csv', skiprows=1)
    mu_p = np.interp(E/1000, xcom[:,0], xcom[:,1]) 
    if mat == 'dy':
        density = 8.55
    elif mat == 'ho':
        density = 8.80
    
    for i in range(N):
        # Generate start position
        start = np.array([side_l*np.random.random() - side_l/2, side_l*np.random.random() - side_l/2, thickness*np.random.random()])
        # Generate a random directional vector
        theta = np.random.uniform(0, 2 * np.pi)  # azimuthal angle
        phi = np.arccos(-np.random.uniform(0, 1))  # polar angle for downward hemisphere
        vec = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])

        # Check whether it will be possible to enter the detector
        det_vec = start + np.abs((det_dist+start[2])/vec[2])*vec # vector to detector plane
        if ((det_vec[0]**2 + det_vec[1]**2) > det_R**2): # miss
            pass
        elif attenuate == False: # Hit if not attenuating
            counts += 1
        else: # Check if attenuated
            bottom_pos = start + np.abs(start[2]/vec[2])*vec
            path = bottom_pos - start
            distance = np.sqrt(path[0]**2 + path[1]**2 + path[2]**2)
            survival_prob = np.exp(-mu_p*density*distance)
            if survival_prob > np.random.random():
                counts+=1

    return counts

def foil_calc_cir(radius, thickness, det_dist, det_R, N, attenuate = False, mat='dy', E=100.10595):

    counts = 0
    paths = np.zeros(shape=(N,3))
    tally = [0,0,0,0,0,0]
    xcom=np.loadtxt(f'xcom_{mat}.csv', skiprows=1)
    mu_p = np.interp(E/1000, xcom[:,0], xcom[:,1]) 
    if mat == 'dy':
        density = 8.55
    elif mat == 'ho':
        density = 8.80
    elif mat == 'mnal':
        density = 2.7
    elif mat == 'aual':
        density = 2.7

    for i in range(N):
        while True:
            # Generate start position
            start = np.array([2*radius*np.random.random() - radius, 2*radius*np.random.random() - radius, thickness*np.random.random()])
            if (start[0]**2 + start[1]**2) <= radius**2:
                break
        # Generate a random directional vector
        theta = np.random.uniform(0, 2 * np.pi)  # azimuthal angle
        phi = np.arccos(-np.random.uniform(0, 1))  # polar angle for downward hemisphere
        vec = np.array([np.sin(phi)*np.cos(theta), np.sin(phi)*np.sin(theta), np.cos(phi)])

        # Check whether it will be possible to enter the detector
        det_vec = start + np.abs((det_dist+start[2])/vec[2])*vec # vector to detector plane
        if ((det_vec[0]**2 + det_vec[1]**2) > det_R**2):
            # print('miss', det_vec)
            pass
        elif attenuate == False:
                counts += 1
        else:
            bottom_pos = start + np.abs(start[2]/vec[2])*vec
            path = bottom_pos - start
            distance = np.sqrt(path[0]**2 + path[1]**2 + path[2]**2)
            survival_prob = np.exp(-mu_p*density*distance)
            if survival_prob > np.random.random():
                counts+=1

    return counts


det_dist = 2
# Ta-182 energies
# energies = [67.6784, 100.10595, 113.6717, 152.42991, 156.3864, 179.39381, 198.35187, 222.1085, 229.3207, 264.074 ,1121.29, 1189.04, 1221.395, 1231.004]
energies = [86.76]
corrections = []
iterations = 1E6
mat = 'ho'
thickness = 0.01
side = 2.5
det_R = 2.465 # 2.335
mnal_radius = 0.6

for energy in energies:

    foil_counts = foil_calc_sq(2.5, thickness, det_dist, 2.465, iterations, mat=mat,  E=energy)
    # foil_counts = foil_calc_cir(mnal_radius, thickness, det_dist, det_R, iterations, mat=mat,  E=energy)
    point_counts = point_calc(det_dist, det_R, iterations)
    corrections.append(foil_counts/point_counts)
    print(f'E = {energy} keV {det_dist} cm {mat} correction factor: {corrections[-1]:.2f}')
    # print(f'Foil: {foil_counts}\nPoint: {point_counts}\nRatio: {foil_counts/point_counts:.3f}')
    # print(f'Point counts: {point_counts} / 10000000 = {point_counts/10E6:.4f}')
    # print(f'Foil counts: {foil_counts} / 10000000 = {foil_counts/10E6:.4f}')

# np.savetxt(f'{mat}_corrections_{det_dist}cm.txt', corrections)