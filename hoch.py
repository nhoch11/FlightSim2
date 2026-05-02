import numpy as np
import json
import matplotlib.pyplot as plt

# Constants
PI = 3.1415926535897932384626433832795
TOLERANCE = 1.0e-14

   

def quat_mult(A, B):

    """
    Multiply two quaternions A and B.
    
    Parameters
    ----------
    A : array-like, shape (4,)
    B : array-like, shape (4,)
    
    Returns
    -------
    quat_result : ndarray, shape (4,)
    """
    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)
    
    quat_result = np.zeros(4, dtype=float)

    quat_result[0] = A[0]*B[0] - A[1]*B[1] - A[2]*B[2] - A[3]*B[3]
    quat_result[1] = A[0]*B[1] + A[1]*B[0] + A[2]*B[3] - A[3]*B[2]
    quat_result[2] = A[0]*B[2] - A[1]*B[3] + A[2]*B[0] + A[3]*B[1]
    quat_result[3] = A[0]*B[3] + A[1]*B[2] - A[2]*B[1] + A[3]*B[0]

    return quat_result



def quat_base_to_dependent(base, quat): 

    base = np.asarray(base, dtype = float)
    quat = np.asarray(quat, dtype = float)
    
    T = np.zeros(4, dtype=float)
    dependent = np.zeros(3, dtype = float)

    T[0] = -base[0]*quat[1] - base[1]*quat[2] - base[2]*quat[3]
    T[1] =  base[0]*quat[0] + base[1]*quat[3] - base[2]*quat[2]
    T[2] = -base[0]*quat[3] + base[1]*quat[0] + base[2]*quat[1]
    T[3] =  base[0]*quat[2] - base[1]*quat[1] + base[2]*quat[0]
    
    dependent[0] = quat[0]*T[1] - quat[1]*T[0] - quat[2]*T[3] + quat[3]*T[2]
    dependent[1] = quat[0]*T[2] + quat[1]*T[3] - quat[2]*T[0] - quat[3]*T[1]
    dependent[2] = quat[0]*T[3] - quat[1]*T[2] + quat[2]*T[1] - quat[3]*T[0]

    return dependent


def quat_dependent_to_base(dependent, quat):

    dependent = np.asarray(dependent, dtype=float)
    quat = np.asarray(quat, dtype=float)
    T = np.zeros(4, dtype=float)
    base = np.zeros(3, dtype=float)
    
    T[0] =  dependent[0]*quat[1] + dependent[1]*quat[2] + dependent[2]*quat[3]
    T[1] =  dependent[0]*quat[0] - dependent[1]*quat[3] + dependent[2]*quat[2]
    T[2] =  dependent[0]*quat[3] + dependent[1]*quat[0] - dependent[2]*quat[1]
    T[3] = -dependent[0]*quat[2] + dependent[1]*quat[1] + dependent[2]*quat[0]
    
    base[0] = quat[0]*T[1] + quat[1]*T[0] + quat[2]*T[3] - quat[3]*T[2]
    base[1] = quat[0]*T[2] - quat[1]*T[3] + quat[2]*T[0] + quat[3]*T[1]
    base[2] = quat[0]*T[3] + quat[1]*T[2] - quat[2]*T[1] + quat[3]*T[0]

    return base


def quat_norm(quat):

    quat = np.asarray(4, dtype=float)

    norm = np.sqrt(sum(quat**2))
    quat = quat/norm

    return quat



# euler angles in radians
def euler_to_quat(euler):

    euler = np.asarray(euler, dtype=float)
    quat = np.zeros(4, dtype=float)

    C_half_phi = np.cos(euler[0]/2.) 
    S_half_phi = np.sin(euler[0]/2.) 
    C_half_theta = np.cos(euler[1]/2.)
    S_half_theta = np.sin(euler[1]/2.)
    C_half_psi = np.cos(euler[2]/2.)
    S_half_psi = np.sin(euler[2]/2.)

    quat[0] = C_half_phi*C_half_theta*C_half_psi + S_half_phi*S_half_theta*S_half_psi
    quat[1] = S_half_phi*C_half_theta*C_half_psi - C_half_phi*S_half_theta*S_half_psi
    quat[2] = C_half_phi*S_half_theta*C_half_psi + S_half_phi*C_half_theta*S_half_psi
    quat[3] = C_half_phi*C_half_theta*S_half_psi - S_half_phi*S_half_theta*C_half_psi
    
    return quat


def quat_to_euler(quat):
    
    quat = np.asarray(quat, dtype=float)
    euler = np.zeros(3, dtype=float)

    a = quat[0]*quat[2] - quat[1]*quat[3]
    
    if (abs(a - 0.5) < 1.e-12):
        euler[0] = 2.0*np.arcsin(quat[1]/np.cos(PI/4.0))
        euler[1] = PI/2.0
        euler[2] = 0.
    elif (abs(a + 0.5) < 1.e-12):
        euler[0] = 2.0*np.arcsin(quat[1]/np.cos(PI/4.0));
        euler[1] = -PI/2.0
        euler[2] = 0.0
    else:
        euler[0] = np.arctan2( 2.0*(quat[0]*quat[1] + quat[2]*quat[3]), quat[0]*quat[0] + quat[3]*quat[3] - quat[1]*quat[1] - quat[2]*quat[2])
        euler[1] = np.arcsin(2.0*(quat[0]*quat[2] - quat[1]*quat[3]))
        euler[2] = np.arctan2( 2.0*(quat[0]*quat[3] + quat[1]*quat[2]), quat[0]*quat[0] + quat[1]*quat[1] - quat[2]*quat[2] - quat[3]*quat[3])
    
    return euler


def cross3(A, B):

    A = np.asarray(A, dtype=float)
    B = np.asarray(B, dtype=float)

    ans = np.zeros(3, dtype=float)

    ans[0] =   A[1]*B[2] - A[2]*B[1]
    ans[1] =   A[2]*B[0] - A[0]*B[2]
    ans[2] =   A[0]*B[1] - A[1]*B[0]

    return ans


def plot_view_plane(y_c_vp, z_c_vp, aspect_ratio_vp, grid_color, y_vp, z_vp):

    fig = plt.figure(figsize=(aspect_ratio_vp*5.0,5.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
    plt.axis("off")
    ax.axes.set_xlim( y_c_vp[0], y_c_vp[2])
    ax.axes.set_ylim(-z_c_vp[1],-z_c_vp[0])
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axes.set_aspect("equal")
    # ax.scatter(y_vp, -z_vp, color = "k", s = 4)
    for i in range(0,len(y_vp), 2):
        ax.plot([y_vp[i], y_vp[i+1]], [-z_vp[i], -z_vp[i+1]], color = grid_color)
    fig.canvas.draw()
    plt.show()


def init_view_plane(d_vp, angle_vp, aspect_ratio_vp):

    # width of view plane
    w_vp = 2*d_vp*np.tan(angle_vp)

    # height of view plane
    h_vp = w_vp/aspect_ratio_vp

    # corner points in camera coordinate system
    x_c_vp = np.zeros(4)
    y_c_vp = np.zeros(4)
    z_c_vp = np.zeros(4)

    x_c_vp[:] = d_vp
    
    y_c_vp[0:2] = -0.5*w_vp 
    y_c_vp[2:4] =  0.5*w_vp 
    
    z_c_vp[[0,3]] = -0.5*h_vp 
    z_c_vp[[1,2]] =  0.5*h_vp 

    return x_c_vp, y_c_vp, z_c_vp

    # print("x_c_vp = ", x_c_vp)
    # print("y_c_vp = ", y_c_vp)
    # print("z_c_vp = ", z_c_vp)


def convert_camera_to_earth_fixed(point_cam, cam_pos_earth, cam_quat):

    point_cam = np.asarray(point_cam, dtype=float) # point in camera coordinates
    cam_pos_earth = np.asarray(cam_pos_earth, dtype = float) # camera position in earth fixed coordinates
    cam_quat = np.asarray(cam_quat, dtype=float) # quat of euler angles of camera relative to earth fixed 

    point_earth = quat_dependent_to_base(point_cam, cam_quat) + cam_pos_earth

    return point_earth


def generate_ground_plane(input_json):

    # read json
    json_string=open(input_json).read()
    json_vals = json.loads(json_string)

    # surface grid
    P_fc = json_vals["camera"]["location[ft]"]
    eul_c = json_vals["camera"]["orientation[deg]"]
    eul_c = np.radians(eul_c)
    d_vp = json_vals["camera"]["view_plane"]["distance[ft]"]
    angle_vp = json_vals["camera"]["view_plane"]["angle[deg]"]
    angle_vp = np.radians(angle_vp)
    aspect_ratio_vp = json_vals["camera"]["view_plane"]["aspect_ratio"]

    alt = json_vals["ground"]["altitude[ft]"]
    grid_number = json_vals["ground"]["grid_number"]
    grid_scale = json_vals["ground"]["grid_scale[ft]"]
    grid_color = json_vals["ground"]["color"]

    grid_points_earth_x = []
    grid_points_earth_y = []
    grid_points_earth_z = []

    for i in reversed(range(grid_number)):
        grid_points_earth_x.append(-grid_number*grid_scale)
        grid_points_earth_x.append(grid_number*grid_scale)
        grid_points_earth_y.append(-(i + 1)*grid_scale)
        grid_points_earth_y.append(-(i + 1)*grid_scale)
        grid_points_earth_z.append(-alt)
        grid_points_earth_z.append(-alt)

    grid_points_earth_x.append(-grid_number*grid_scale)
    # grid_points_earth_x.append(-grid_number*grid_scale)
    grid_points_earth_y.append(0)
    # grid_points_earth_y.append(0)
    grid_points_earth_z.append(-alt)
    # grid_points_earth_z.append(-alt)

    # do the first half of the lists
    grid_points_earth_x.extend([-x for x in grid_points_earth_x])
    grid_points_earth_y.extend([-y for y in reversed(grid_points_earth_y)])
    grid_points_earth_z.extend(grid_points_earth_z)
    x_copy = grid_points_earth_x.copy()

    # then append x to y and y to x
    grid_points_earth_x.extend(grid_points_earth_y)
    grid_points_earth_y.extend(x_copy)
    grid_points_earth_z.extend(grid_points_earth_z)

    # print("x = ", grid_points_earth_x)
    # print("y = ", grid_points_earth_y)
    # print("z = ", grid_points_earth_z)

    # init_view_plane
    x_c_vp, y_c_vp, z_c_vp = init_view_plane(d_vp, angle_vp, aspect_ratio_vp)

    # convert to eath fixed
    x_f_vp = np.zeros(4)
    y_f_vp = np.zeros(4)
    z_f_vp = np.zeros(4) 

    quat = euler_to_quat(eul_c)

    for i in range(4):
        [x_f_vp[i], y_f_vp[i], z_f_vp[i]] = convert_camera_to_earth_fixed([x_c_vp[i], y_c_vp[i], z_c_vp[i]], P_fc, quat)
            

    # calc P0
    P0 = np.asarray([np.average(x_f_vp), np.average(y_f_vp), np.average(z_f_vp)])
    print("P0 = ", P0)

    # calc n_f_vp
    n_f_vp = cross3([x_f_vp[0] - P0[0], y_f_vp[0] - P0[1], z_f_vp[0] - P0[2]], [x_f_vp[1] - P0[0], y_f_vp[1] - P0[1], z_f_vp[1] - P0[2]])
    print("n_f_vp = ", n_f_vp)

    print("P0 - Pc = ", P0-P_fc)

    l_cax = np.asarray(grid_points_earth_x) - P_fc[0] 
    l_cay = np.asarray(grid_points_earth_y) - P_fc[1]
    l_caz = np.asarray(grid_points_earth_z) - P_fc[2]

    l_ca = np.stack((l_cax, l_cay, l_caz), axis=1)

    print("l_cax = ", l_cax)
    print("l_cay = ", l_cay)
    print("l_caz = ", l_caz)

    gamma = np.dot((P0 - P_fc), n_f_vp)/np.dot(l_ca, n_f_vp)
    print("gamma = ", gamma)


    vp_points = np.zeros((len(l_ca),3))

    for i in range(len(l_ca)):
        vp_points[i] = quat_base_to_dependent(gamma[i]*l_ca[i,:], quat)

    print("x in plot = ", vp_points[:,1])
    print("y in plot = ", -vp_points[:,2])


    plot_view_plane(y_c_vp, z_c_vp, aspect_ratio_vp, grid_color, vp_points[:,1], vp_points[:,2])

    

