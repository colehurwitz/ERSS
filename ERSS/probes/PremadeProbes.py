from BaseProbe import BaseProbe
import numpy as np
np.set_printoptions(suppress=True)
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d.art3d import Poly3DCollection, Line3DCollection
from mpl_toolkits.mplot3d import Axes3D
from collections import namedtuple
import pandas as pd
from scipy.spatial import cKDTree
import math
Channel = namedtuple("Channel", ["id", "pos", "type"])

class NeuropixelsProbe(BaseProbe):
    '''Dense probe that can record with up to 384 (10 reference) channels. 1 cm shank length and 70 micron width.


    Attributes:
        width (float)                   Width of the probe in microns
        num_channels (int)              Number of channels on the probe
        origin (numpy.ndarray (3D))     The point where the probe starts in the 3D environment
        reference_channels (array_like) channel ids for all reference channels on the given channels
        channels_df (DataFrame)         The dataframe containing all channels on the probe
    
    '''
    def __init__(self, width=70.0, num_channels=384, origin=None, reference_channels=None, channel_df=None):
        BaseProbe.__init__(self, num_channels, origin, channel_df)
        
        if(self.num_channels > 384):
            raise ValueError('Number of recording channels cannot be larger than 384 for a Neuropixels probe')
        
        self.width = width         
        if(reference_channels is not None):
            self.reference_channels = reference_channels
        else:
            reference_channels = np.random.choice(self.num_channels, int(10*(self.num_channels/384.0)), replace=False)
            self.reference_channels = reference_channels
        if(origin is None):
            self.origin = np.asarray([0, 0, 0])
        else:
            self.origin = origin
        depth  = self.num_channels*10 + 20
        vertices = np.array([[-self.width/2, 0, 0], [self.width/2, 0, 0], [-self.width/2, 0, depth],  
                         [self.width/2, 0, depth]]) + self.origin
        self.vertices = vertices
        channels = []
        initial_channels_pos = np.asarray([[-self.width/2 + 11, 0, 20], [-self.width/2 + 43, 0, 20], 
                                           [-self.width/2 + 27, 0, 40], [-self.width/2 + 59, 0, 40]]) + self.origin
        more_channels = True
        z_offset = 0
        i = 0
        while more_channels:
            for channel_pos in initial_channels_pos:
                channel_pos = channel_pos + np.asarray([0,0, z_offset])
                if(i in self.reference_channels):
                    channel = Channel(i, channel_pos, 'ref')
                    i += 1
                else:
                    channel = Channel(i, channel_pos, 'rec')
                    i += 1
                channels.append(channel)
                if len(channels) == self.num_channels:
                    more_channels = False
                    break
            z_offset += 40       
        channel_df = pd.DataFrame(channels, columns=Channel._fields)
        self.channel_df = channel_df
        
    def draw(self, ax=None, excluded_channel_ids=None):
        if(ax is None):
            fig = plt.figure(figsize=(30,30))
            ax = fig.add_subplot(111, projection='3d')
            
        sides = [[self.vertices[0],self.vertices[2]], [self.vertices[0],self.vertices[1]], 
                 [self.vertices[2],self.vertices[3]], [self.vertices[1],self.vertices[3]]]
        ax.add_collection3d(Poly3DCollection(sides, facecolors='red', linewidths=1, edgecolors='r'))
        
        if(excluded_channel_ids is None):
            for i, row in self.channel_df.iterrows():
                if(row['type'] == 'ref'):
                    channel = Channel(i, row['pos'], 'ref')
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='orange')
                else:
                    channel = Channel(i, row['pos'], 'rec')
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='black')
        else:
            for i, row in self.channel_df.iterrows():
                if(i not in excluded_channel_ids):
                    if(row['type'] == 'ref'):
                        channel = Channel(i, row['pos'], 'ref')
                        ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='orange')
                    else:
                        channel = Channel(i, row['pos'], 'rec')
                        ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='black')
                else:
                    continue
        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        ax.set_xlabel('x (microns)')
        ax.set_ylabel('y (microns)')
        ax.set_zlabel('Depth (Microns)')

    def rotate(self, theta, axis, ax=None, plot=False):
        if(ax is None and plot):
            fig = plt.figure(figsize=(30,30))
            ax = fig.add_subplot(111, projection='3d')
            
        for i, vertice in enumerate(self.vertices):
            self.vertices[i] = np.dot(rotation_matrix(axis, theta), vertice)
        sides = [[self.vertices[0],self.vertices[2]], [self.vertices[0],self.vertices[1]], 
                 [self.vertices[2],self.vertices[3]], [self.vertices[1],self.vertices[3]]]
        if(plot):
            ax.add_collection3d(Poly3DCollection(sides, facecolors='red', linewidths=1, edgecolors='r'))
        
        new_channels = []
        for i, row in self.channel_df.iterrows():
            if(row['type'] == 'ref'):
                channel = Channel(i, np.dot(rotation_matrix(axis, theta), row['pos']), 'ref')
                if(plot):
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='orange')
                new_channels.append(channel)
            else:
                channel = Channel(i, np.dot(rotation_matrix(axis, theta), row['pos']), 'rec')
                if(plot):
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='black')
                new_channels.append(channel)
                
        channel_df = pd.DataFrame(new_channels, columns=Channel._fields)
        self.channel_df = channel_df
        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        ax.set_xlabel('x (microns)')
        ax.set_ylabel('y (microns)')
       
    def shift(self, dist, axis, ax, plot):
        if(ax is None and plot):
            fig = plt.figure(figsize=(30,30))
            ax = fig.add_subplot(111, projection='3d')
            
        for i, vertice in enumerate(self.vertices):
            self.vertices[i] = vertice + axis*dist
        sides = [[self.vertices[0],self.vertices[2]], [self.vertices[0],self.vertices[1]], 
                 [self.vertices[2],self.vertices[3]], [self.vertices[1],self.vertices[3]]]
        if(plot):
            ax.add_collection3d(Poly3DCollection(sides, facecolors='red', linewidths=1, edgecolors='r'))
        
        new_channels = []
        for i, row in self.channel_df.iterrows():
            if(row['type'] == 'ref'):
                channel = Channel(i, row['pos'] + axis*dist, 'ref')
                if(plot):
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='orange')
                new_channels.append(channel)
            else:
                channel = Channel(i, row['pos'] + axis*dist, 'rec')
                if(plot):
                    ax.scatter3D(channel.pos[0], channel.pos[1], channel.pos[2], c='black')
                new_channels.append(channel)
                
        channel_df = pd.DataFrame(new_channels, columns=Channel._fields)
        self.channel_df = channel_df
        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        ax.set_xlabel('x (microns)')
        ax.set_ylabel('y (microns)')
        ax.set_zlabel('Depth (Microns)')    
        
    
def rotation_matrix(axis, theta):
    """
    Return the rotation matrix associated with counterclockwise rotation about
    the given axis by theta radians.
    """
    axis = np.asarray(axis)
    axis = axis / math.sqrt(np.dot(axis, axis))
    a = math.cos(theta / 2.0)
    b, c, d = -axis * math.sin(theta / 2.0)
    aa, bb, cc, dd = a * a, b * b, c * c, d * d
    bc, ad, ac, ab, bd, cd = b * c, a * d, a * c, a * b, b * d, c * d
    return np.array([[aa + bb - cc - dd, 2 * (bc + ad), 2 * (bd - ac)],
                     [2 * (bc - ad), aa + cc - bb - dd, 2 * (cd + ab)],
                     [2 * (bd + ac), 2 * (cd - ab), aa + dd - bb - cc]])
