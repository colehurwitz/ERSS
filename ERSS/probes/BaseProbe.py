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

class BaseProbe(object):
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
