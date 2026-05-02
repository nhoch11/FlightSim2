import numpy as np
import json
import sys
import time
import matplotlib.pyplot as plt
import hoch
from connection_m import connection

class Camera:

    def __init__(self, json_data_camera):
        self.loc_f = json_data_camera["location[ft]"]
        eul_c = json_data_camera["orientation[deg]"]
        self.eul_c = np.radians(eul_c)
        self.quat = hoch.euler_to_quat(self.eul_c)
        self.d_vp = json_data_camera["view_plane"]["distance[ft]"]
        self.angle_vp = json_data_camera["view_plane"]["angle[deg]"]
        self.angle_vp = np.radians(self.angle_vp)
        self.aspect_ratio_vp = json_data_camera["view_plane"]["aspect_ratio"]

        # initi view plane
        # width of view plane
        self.w_vp = 2*self.d_vp*np.tan(self.angle_vp)

        # height of view plane
        self.h_vp = self.w_vp/self.aspect_ratio_vp

        # corner points in camera coordinate system
        self.x_vp_cam = np.zeros(4)
        self.y_vp_cam = np.zeros(4)
        self.z_vp_cam = np.zeros(4)

        self.x_vp_cam[:] = self.d_vp
        
        self.y_vp_cam[0:2] = -0.5*self.w_vp 
        self.y_vp_cam[2:4] =  0.5*self.w_vp 
        
        self.z_vp_cam[[0,3]] = -0.5*self.h_vp 
        self.z_vp_cam[[1,2]] =  0.5*self.h_vp 

        # print("x_vp_cam = ", self.x_vp_cam)
        # print("y_vp_cam = ", self.y_vp_cam)
        # print("z_vp_cam = ", self.z_vp_cam)

        # set up earth fixed vp loc
        self.x_vp_f = np.zeros(4)
        self.y_vp_f = np.zeros(4)
        self.z_vp_f = np.zeros(4) 

        self.lines2D = np.zeros((2,2))

    def set_state(self, location, quat):

        self.loc_f = location
        self.quat = quat

        # convert from camera coord to earth fixed
        for i in range(4):
            [self.x_vp_f[i], self.y_vp_f[i], self.z_vp_f[i]] = hoch.quat_dependent_to_base([self.x_vp_cam[i], self.y_vp_cam[i], self.z_vp_cam[i]], self.quat) + self.loc_f

        # print("\nx_vp_f = ", self.x_vp_f)
        # print("y_vp_f = ", self.y_vp_f)
        # print("z_vp_f = ", self.z_vp_f)

        # calc P0
        self.P0 = np.asarray([np.average(self.x_vp_f), np.average(self.y_vp_f), np.average(self.z_vp_f)])
        # print("P0 = ", self.P0)

        # calc n_f_vp
        self.n_vp_f = np.cross([self.x_vp_f[0] - self.P0[0], self.y_vp_f[0] - self.P0[1], self.z_vp_f[0] - self.P0[2]], [self.x_vp_f[1] - self.P0[0], self.y_vp_f[1] - self.P0[1], self.z_vp_f[1] - self.P0[2]])
        # print("n_vp_f = ", self.n_vp_f)

        # print("P0 - Pc = ", self.P0-self.loc_f)


class LinesObject:

    def __init__(self, json_data, ax):

        if json_data["type"] == "grid":

            ground_alt = json_data["altitude[ft]"]
            grid_number = json_data["grid_number"]
            grid_scale = json_data["grid_scale[ft]"]
            grid_color = json_data["color"]

            num_x_lines = grid_number*2 + 1
            self.num_lines = num_x_lines*2
            self.num_points = 2*self.num_lines

            self.points = np.zeros((2*self.num_lines,3))
            self.lines = np.zeros((self.num_lines,2), dtype=int)

            for i in range(num_x_lines):
                # lines parallel with x axis, moving from -y to y
                self.points[2*i,:]     = [-grid_number*grid_scale, (i - grid_number)*grid_scale, -ground_alt] # left point
                self.points[2*i + 1,:] = [ grid_number*grid_scale, (i - grid_number)*grid_scale, -ground_alt] # right point

                # store indices of the points that make up a line
                self.lines[i,0] = 2*i
                self.lines[i,1] = 2*i + 1

            for i in range(num_x_lines):
                # lines parallel with y axis, moving from -x to x
                self.points[2*i + self.num_lines,:]     = [(i - grid_number)*grid_scale, -grid_number*grid_scale, -ground_alt] # left point
                self.points[2*i + self.num_lines + 1,:] = [(i - grid_number)*grid_scale,  grid_number*grid_scale, -ground_alt] # right point

                # store indices of the points that make up a line
                self.lines[i + num_x_lines,0] = 2*i + self.num_lines
                self.lines[i + num_x_lines,1] = 2*i + self.num_lines + 1

            # print("points = ", self.points)


            # set up ax plot
            self.ax, = ax.plot([], [], linestyle = "-", color = grid_color)

            self.points2D = np.zeros((self.num_points, 2))
            self.l_ca = np.zeros((self.num_points,3))
            self.gamma = np.zeros(self.num_points)
            self.pf = np.zeros((self.num_points, 3))
            
            self.lines2D = np.full((self.num_lines*3, 2), None, dtype = object)
            # self.lines2D = np.full((2 * self.num_lines, 2), np.nan, dtype=float)

    def draw(self, camera):

        self.l_ca[:] = self.points[:] - camera.loc_f
        # print("l_ca = ", self.l_ca)

        

        gamma = np.dot((camera.P0 - camera.loc_f), camera.n_vp_f)/np.dot(self.l_ca, camera.n_vp_f)
        # print("gamma = ", gamma)

        self.pf[:] = self.l_ca[:]*gamma[:,None]


        vp_points = np.zeros((len(gamma),3))

        for i in range(len(gamma)):
            vp_points[i] = hoch.quat_base_to_dependent(gamma[i]*self.l_ca[i,:], camera.quat)
        
        self.points2D[:,0] =  vp_points[:,1]
        self.points2D[:,1] = -vp_points[:,2]

        # print("x = ", self.points2D[:,0])
        # print("y = ", self.points2D[:,1])
        zer = np.zeros(2)
        # print("check 123")
        for i in range(self.num_lines):
            i0 = self.lines[i,0]
            i1 = self.lines[i,1]
            
            mid = 0.5 * (self.points[i0] + self.points[i1])
            d_mid = np.linalg.norm(mid - camera.loc_f)

            if d_mid > 30000:
                self.lines2D[3*i, :] = zer
                self.lines2D[3*i + 1, :] = zer
                continue

            
            if (gamma[i0]>0. and gamma[i1]>0.): # if entire line is in front of the view plane
                # print("check0")
                self.lines2D[3*i,:]   = self.points2D[i0,:]
                self.lines2D[3*i+1,:] = self.points2D[i1,:]
                # print(self.lines2D[3*i,:])
                # print(self.lines2D[3*i+1,:])
            elif (gamma[i0]<0. and gamma[i1]<0): # entire line is behind the viewer
                # print("check1")
                self.lines2D[3*i,:]   = zer
                self.lines2D[3*i+1,:] = zer
                # print(self.lines2D[3*i,:])
                # print(self.lines2D[3*i+1,:])

            else:
                # print("check2")
                # print(gamma[i0])
                # print(gamma[i1])

                l = self.points[i1] - self.points[i0]
                gam = np.dot(camera.P0 - self.points[i0], camera.n_vp_f)/ np.dot(l,camera.n_vp_f)
                clipped_3d = hoch.quat_base_to_dependent(self.points[i0] + gam*l - camera.loc_f, camera.quat)
                clipped_2d = np.asarray([clipped_3d[1], -clipped_3d[2]], dtype=float)
            

                if (gamma[i0]<0.):
                    # print("check3")
                    self.lines2D[3*i,:]   = clipped_2d
                    self.lines2D[3*i+1,:] = self.points2D[i1,:]

                    # print(self.lines2D[3*i,:])
                    # print(self.lines2D[3*i+1,:])

                if (gamma[i1]<0.):
                    # print("check5")
                    # l = self.points[i0] - self.points[i1]
                    # gam = np.dot(camera.P0 - self.points[i1], camera.n_vp_f)/ np.dot(l,camera.n_vp_f)
                    # clipped_3d = hoch.quat_base_to_dependent(gam*l, camera.quat)
                    # clipped_2d = np.asarray([clipped_3d[1], -clipped_3d[2]], dtype=float)
                    self.lines2D[3*i,:]   = self.points2D[i0,:]
                    self.lines2D[3*i+1,:] = clipped_2d

                    # print(self.lines2D[3*i,:])
                    # print(self.lines2D[3*i+1,:])
                # print("neither")
        
        # print(" x = ", self.lines2D[:,0])
        # print(" y = ", self.lines2D[:,1])
        self.ax.set_data(self.lines2D[:,0], self.lines2D[:,1])

    


if __name__ == "__main__": 

    input_json = "view_input.json"

    # read json
    json_string=open(input_json).read()
    json_data = json.loads(json_string)

    # create connection
    states_connection = connection(json_data["connections"]["receive_states"])

    # init camera
    camera = Camera(json_data["camera"])

    # init plot
    fig = plt.figure(figsize=(camera.aspect_ratio_vp*5.0,5.0))
    ax = fig.add_subplot(111)
    plt.subplots_adjust(top=1.0, bottom=0.0, left=0.0, right=1.0)
    plt.axis("off")
    ax.axes.set_xlim( camera.y_vp_cam[0], camera.y_vp_cam[2])
    ax.axes.set_ylim(-camera.z_vp_cam[1],-camera.z_vp_cam[0])
    ax.axes.xaxis.set_ticklabels([])
    ax.axes.yaxis.set_ticklabels([])
    ax.set_xticks([])
    ax.set_yticks([])
    ax.axes.set_aspect("equal")
    fig.canvas.draw()
    plt.show(block=False)

    # initialize ground plane object
    ground = LinesObject(json_data["scene"]["ground"], ax)

    # initalize location
    location = camera.loc_f
    quat = camera.quat

    fps = 0
    frame = 0

    # while(location[0]<0.0):
    #     time_begin = time.time()
    #     camera.set_state(location, quat)
    #     ground.draw(camera)
    #     fig.canvas.draw()
    #     fig.canvas.flush_events()
    #     location[0] += 0.1
    #     time_end = time.time()
    #     fps = 1/(time_end-time_begin)
    #     print("graphics rate = ", fps)
        # plt.show(block=True)
        # sys.exit()

    while(frame < 1000):
        time_begin = time.time()

        # get updated state from physics
        states = states_connection.recv()
        # print("states = ",states)
        location = states[7:10]
        # print("location = ", location)
        quat = states[10:]
        # print("quat = ", quat)
        camera.set_state(location, quat)
        # print("after set state")
        ground.draw(camera)
        fig.canvas.draw()
        fig.canvas.flush_events()
        time_end = time.time()
        # frame += 1
        fps = 1/(time_end-time_begin)
        print("graphics rate = ", fps)