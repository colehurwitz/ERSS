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

class NeuralProbe(object):
    '''A probe that is used for recording extracellular activity.

    Attributes:
        num_channels (int)           Number of channels on the probe
        origin (numpy.ndarray (3D))  The point where the probe starts in the 3D environment
        channels_df (DataFrame)      The dataframe containing all channels on the probe
        
    '''
    def __init__(self, num_channels, origin, channels_df):
        self.num_channels = num_channels
        self.origin = origin
        self.channel_df = channels_df
        
    def draw(self, ax, excluded_channel_ids=None):
        '''Draws the probe object in 3D space.

        Parameters
        ----------
        ax: Axes3D
            The 3D axes for plotting the probe
        excluded_channel_ids: array_like
            list of channel ids to be excluded from the probe during the draw step

        '''
        raise NotImplementedError("The draw function is not implemented for \
            this probe")
    
    def rotate(self, theta, axis, ax, plot):
        '''Rotates the probe object around an axis with a given angle. Plotting is optional,

        Parameters
        ----------
        theta: float
            The amount of radians to rotate the probe.
        axis: numpy.ndarray (3D)
            The 3D vector the probe is rotated around
        ax: Axes3D
            The 3D axes for plotting the probe

        '''
        raise NotImplementedError("The rotate function is not implemented for \
            this probe")
    
    def shift(self, dist, axis, ax, plot):
        '''Shifts the probe object in the direction of the axis a given dist. Plotting is optional,

        Parameters
        ----------
        dist: float
            The amount of microns to shift the probe
        axis: numpy.ndarray (3D)
            The 3D vector direction in which the probe is shifted
        ax: Axes3D
            The 3D axes for plotting the probe

        '''
        raise NotImplementedError("The shift function is not implemented for \
            this probe")

    def getChannels(self, excluded_channel_ids=None):
        '''Returns DataFrame of all channels in the probe

        Returns
        -------
        self.channels_df: DataFrame
            DataFrame containing all channels in the probe
        '''
        all_channel_ids_set = set(range(self.num_channels))
        excluded_channel_ids_set = set(excluded_channel_ids)
        recording_channels_indices_set = all_channel_ids_set - excluded_channel_ids_set
        recording_channels_indices = list(recording_channels_indices_set)
        recording_channels_df = self.channel_df.loc[recording_channels_indices].copy()
        
        return recording_channels_df
    
    def getNeuronsRadius(self, neuron_df, radius=100.0, excluded_channel_ids=None):
        '''Returns DataFrame of all neurons in a radius around the channels of the probe (sans excluded channels)

        Parameters
        ----------
        neuron_df: DataFrame
            DataFrame containing all neurons in the brain region
        radius: float
            Radius around channels in which neurons are returned
        excluded_channel_ids: array_like
            list of channel ids to be excluded from the probe during the getting neurons step

        Returns
        -------
        close_neuron_df: DataFrame
            DataFrame containing all neurons within the given radius of any channel
        '''
        close_neurons_indices = set()
        
        #Create kdtree of neuron positions for fast look up
        neuron_pos = np.asarray(list(neuron_df['pos']))
        tree = cKDTree(neuron_pos)
        for i, row in self.channel_df.iterrows():
            if(excluded_channel_ids is None):
                channel_pos = row['pos']
                indices = tree.query_ball_point(channel_pos, radius)        
                for idx in indices:
                    close_neurons_indices.add(idx)
            else:
                if(i not in excluded_channel_ids):
                    channel_pos = row['pos']
                    indices = tree.query_ball_point(channel_pos, radius)        
                    for idx in indices:
                        close_neurons_indices.add(idx)
                else:
                    continue
        close_neurons_indices = list(close_neurons_indices)
        close_neuron_df = neuron_df.loc[close_neurons_indices].copy()
        
        return close_neuron_df


class NeuroPixel(NeuralProbe):
    '''Dense probe that can record with up to 384 (10 reference) channels. 1 cm shank length and 70 micron width.


    Attributes:
        width (float)                   Width of the probe in microns
        num_channels (int)              Number of channels on the probe
        origin (numpy.ndarray (3D))     The point where the probe starts in the 3D environment
        reference_channels (array_like) channel ids for all reference channels on the given channels
        channels_df (DataFrame)         The dataframe containing all channels on the probe
    
    '''
    def __init__(self, width=70.0, num_channels=384, origin=None, reference_channels=None, channel_df=None):
        NeuralProbe.__init__(self, num_channels, origin, channel_df)
        
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
        ax.set_zlabel('Depth (Microns)')
    
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