import json
from os.path import join
import numpy as np
import socket
import struct
import datetime
#  from datetime import datetime as dt
#  from datetime import timedelta as td

def csvLineWrite(*obj, sep=',', end='\n'):
    line = ''
    k = len(obj)
    if k == 0: return end
    for i,o in enumerate(obj):
        Sep = sep
        if i == k-1: Sep = ''
        if type(o) == float or type(o) == int or type(o) == np.float64:
            line += '{:23.16e}'.format(o) + Sep
        elif o == None:
            line += Sep
        else:
            line += str(o) + Sep
    return line+end

##======================================================================

class database():
    def __init__(self, fn, pn):
        if pn == None:
            input_file = fn
        else:
            input_file = join(pn, fn)
        
        f = open(input_file, 'r')
        data_raw = f.readlines()
        f.close()
        
        self.params = {}
        
        n_rows = len(data_raw)
        n_comments = 0
        n_parameters = 0
        n_pts = 0
        flag = False
        for i, line in enumerate(data_raw):
            
            if line[0] == '#':
                n_comments += 1
                continue
            
            if flag:
                n_pts += 1
                #  cols = csvLineRead(line)
                cols = line[:-len('\n')].split(',')
                if len(cols) != 1+self.n_dv:
                    print(f'Error reading database csv {fn}! Row {i+1} has a different number of columns! Quitting...')
                    exit(1)
                continue
            
            #  cols = csvLineRead(line)
            cols = line[:-len('\n')].split(',')
            
            if cols[0] == 'paramater':
                n_parameters += 1
                self.params[cols[1]] = cols[2]
            else:
                flag = True
                headings_row = i
                self.n_dv = len(cols) - 1
                self.dep_vars = cols[1:]
        
        if n_rows != n_comments + n_parameters + 1 + n_pts:
            print(f'Error reading database csv {fn}! Number of rows, comments, parameters and data points do not agree! Quitting...')
            exit(1)
        
        self.n_pts = n_pts
        self.x = [None] * n_pts
        self.y = np.zeros((n_pts, self.n_dv))
        
        row = -1
        for line in data_raw[headings_row+2:]:
            if line[0] == '#': continue
            row += 1
            cols = line[:-len('\n')].split(',')
            vals = [float(i) for i in cols]
            self.x[row] = vals[0]
            self.y[row,:] = vals[1:]
    
    def interpolate(self, t):
        imax = None
        imin = None
        if t in self.x:
            i = self.x.index(t)
            if i == self.n_pts - 1:
                imax = i
                imin = i-1
            else:
                imin = i
                imax = i+1
        else:
            for i in range(self.n_pts-1):
                if self.x[i] < t and t < self.x[i+1]:
                    imin = i
                    imax = i+1
                    break
        if imax == None: raise ValueError(f'Error interpolating database. Value {t} is not in the bound: ({self.x[0]}, {self.x[-1]}). Quitting...')
        return (self.y[imax,:] - self.y[imin,:]) / (self.x[imax] - self.x[imin]) * (t - self.x[imin]) + self.y[imin,:]

##======================================================================

class channel():
    def __general_init__(self, json_data):
        self.n = json_data['number_of_values']
        temp = json_data['type']
        if temp == 'send':
            self.only_send = True
        elif temp == 'receive':
            self.only_send = False
            self.vals = [0]*self.n
        else:
            print('Error initializing channel. Type {} not recognized! Should be either send or receive. Quitting...')
            exit(1)
        self.verbose = json_data.get('verbose', True)

class file_channel(channel):
    def __init__(self, json_data):
        self.__general_init__(json_data)
        
        temp = json_data['channel_type']
        if temp != 'file':
            print('Error initializing channel. Initializing as file channel but channel_type does not match! Quitting...')
            exit(1)
        
        self.fn = json_data['filename']
        pn = json_data.get('pathname', None)
        
        if self.only_send:
            if pn != None:
                input_file = join(pn, self.fn)
            else:
                input_file = self.fn
            self.fid = open(input_file, 'w')
            cols = json_data.get('column_names', None)
            if cols != None:
                self.fid.write(csvLineWrite(*cols))
        else:
            #  temp = json_data['database_type']
            self.db = database(self.fn, pn)
            if self.n != self.db.n_dv:
                print('Error initializing file receive channel. Database output does not match the expected amount of values! Quitting...')
                exit(1)
    
    def send(self, vals):
        if not self.only_send:
            print(f'Error with file channel to file {self.fn}. Attempting to send data, but channel has been setup as receive only! Quitting...')
            exit(1)
        
        self.fid.write(csvLineWrite(*vals))
        if self.verbose:
            print('Sent values to file {}'.format(self.fn))
            print(vals)
    
    def recv(self, x):
        if self.only_send:
            print(f'Error with file channel from file {self.fn}. Attempting to recv data, but channel has been setup as send only! Quitting...')
            exit(1)
        
        self.vals = self.db.interpolate(x)
        if self.verbose:
            print('Recieved values from file {}'.format(self.fn))
            print(self.vals)
        return self.vals

class udp_channel(channel):
    def __init__(self, json_data):
        self.__general_init__(json_data)
        
        temp = json_data['channel_type']
        if temp != 'udp':
            print('Error initializing channel. Initializing as udp channel but channel_type does not match! Quitting...')
            exit(1)
        
        self.port = json_data['port_ID']
        self.double_precision = json_data.get('double_precision', False)
        self.little_endian = json_data.get('little_endian', True)
        
        self.sock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        
        if self.only_send:
            self.ip = json_data.get('IP_address', '127.0.0.1')
        else:
            if self.double_precision:
                self.__buffer__ = 8 * self.n
            else:
                self.__buffer__ = 4 * self.n
            self.ip = json_data.get('IP_address', '0.0.0.0')
            self.blocking = json_data.get('wait_for_data', False)
            self.sock.bind((self.ip, self.port))
            self.sock.setblocking(self.blocking)
    
    def send(self, vals):
        if not self.only_send:
            print(f'Error with udp channel with port {self.port}. Attempting to send data, but channel has been setup as receive only! Quitting...')
            exit(1)
        if len(vals) != self.n:
            print(f'Error with udp channel with port {self.port}. Attempting to send {len(vals)} data values but excpected {self.n} values! Quitting...')
            exit(1)
        
        if self.double_precision:
            data = np.array(vals, dtype=np.float64)
            if self.little_endian:
                packed_data = struct.pack('<' + 'd' * self.n, *data)
            else:
                packed_data = struct.pack(f"!{self.n}d", *data)
        else:
            data = np.array(vals, dtype=np.float32)
            if self.little_endian:
                packed_data = struct.pack('<' + 'f' * self.n, *data)
            else:
                packed_data = struct.pack(f"!{self.n}f", *data)
        self.sock.sendto(packed_data, (self.ip, self.port))
        if self.verbose:
            print('Sent values to udp port {}'.format(self.port))
            print(vals)
    
    def recv(self):
        if self.only_send:
            print(f'Error with udp channel on port {self.port}. Attempting to recv data, but channel has been setup as send only! Quitting...')
            exit(1)
        
        if self.blocking:
            data, addr = self.sock.recvfrom(self.__buffer__)
            self.sock.setblocking(False)
            try:
                while True:
                    d, a = self.sock.recvfrom(self.__buffer__)
                    data, addr = d, a
            except BlockingIOError:
                pass
            self.sock.setblocking(True)
            if self.double_precision:
                if self.little_endian:
                    floats = struct.unpack('<'+'d'*self.n, data)
                else:
                    floats = struct.unpack(f"!{self.n}d", data)
            else:
                if self.little_endian:
                    floats = struct.unpack('<'+'f'*self.n, data)
                else:
                    floats = struct.unpack(f"!{self.n}f", data)
            self.vals = floats
        else:
            try:
                if self.double_precision:
                    if self.little_endian:
                        while True:
                            data, addr = self.sock.recvfrom(self.__buffer__)
                            floats = struct.unpack('<'+'d'*self.n, data)
                            self.vals = floats
                    else:
                        while True:
                            data, addr = self.sock.recvfrom(self.__buffer__)
                            floats = struct.unpack(f"!{self.n}d", data)
                            self.vals = floats
                else:
                    if self.little_endian:
                        while True:
                            data, addr = self.sock.recvfrom(self.__buffer__)
                            floats = struct.unpack('<'+'f'*self.n, data)
                            self.vals = floats
                    else:
                        while True:
                            data, addr = self.sock.recvfrom(self.__buffer__)
                            floats = struct.unpack(f"!{self.n}f", data)
                            self.vals = floats
            except BlockingIOError:
                pass
        
        if self.verbose:
            print('Recieved values to udp port {}'.format(self.port))
            print(self.vals)
        return self.vals

def new_channel(json_data):
    
    temp = json_data['channel_type']
    if temp == 'udp':
        obj = udp_channel(json_data)
    elif temp == 'file':
        obj = file_channel(json_data)
    else:
        print(f'Error initializing channel. Channel_type {temp} not recognized! Quitting...')
        exit(1)
    
    return obj

##======================================================================

class connection():
    def __init__(self, json_data):
        
        self.ch = new_channel(json_data)
        
        refresh_rate = json_data.get('refresh_rate', 0.)
        if refresh_rate <= 0:
            #  self.refresh_time = 1e-12
            self.refresh_time = datetime.timedelta(seconds=1e-12)
        else:
            self.refresh_time = datetime.timedelta(seconds=1 / refresh_rate)
        
        self.time = datetime.datetime.now()
        self.time -= self.refresh_time          ## guarantees first attempt will trigger an update
    
    def check_refresh(self):
        return datetime.datetime.now() - self.time >= self.refresh_time
    
    def send(self, vals):
        if self.check_refresh():
            self.ch.send(vals)
            self.time += self.refresh_time
    
    def recv(self, x=None):
        if self.check_refresh():
            if x == None:
                vals = self.ch.recv()
            else:
                vals = self.ch.recv(x)
            self.time += self.refresh_time
        else:
            vals = self.ch.vals
        return vals



