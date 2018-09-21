# -*- coding: utf-8 -*-
import numpy as np
np.set_printoptions(suppress=True)
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from scipy import stats
from mpl_toolkits.mplot3d import Axes3D
import mpl_toolkits.mplot3d.art3d as art3d
from matplotlib.patches import Circle
import math
import pandas as pd
from scipy.spatial.distance import squareform, pdist
from scipy.spatial import cKDTree
from collections import namedtuple
import matplotlib.ticker as ticker

Layer = namedtuple("Layer", ["id", "region", "height", "radius", "volume", "neuron_count", "density", "in_ratio", 
                             "layer_density_distributions"])
Neuron = namedtuple("Neuron", ["id", "pos", "type", "layer"])

class Column(object):
    """A column from a cortex (details depend on the dictionary being passed in during initialization)

    Attributes:
        layer_df (array_like)     DataFrame containing all layers in the column                                         
        neuron_df (DataFrame)     Dataframe containing neurons in the column
    """
    def __init__(self, cortex_data_dict, layer_df=None, neuron_df=None, radius=200, min_neuron_dist=15):
        self.layer_df = layer_df
        self.neuron_df = neuron_df
        self.min_neuron_dist = min_neuron_dist
        self.radius = radius
        
        #Private function call to build the column from the cortex column dict
        self.__buildColumn(cortex_data_dict)
        
    def plotLayer(self, layer_id, resolution=100, color='r', alpha=.1, x_center=0.0, y_center=0.0, ax=None):
        '''This plots the given layer and the neurons contained in the layer.

        Parameters
        ----------
        layer_id: int
            The id of the layer to be plotted      
        resolution: int
            Resolution of the layer shell
        color: string
            Color of the layer shell
        alpha: float
            Transparency of the layer shell
        x_center: float
            x coordinate of the center of the layer shell
        y_center: float
            y coordinate of the center of the layer shell
        ax: Axes3D
            The 3D axes for plotting the layer
            
        '''
        depth_start = self.depths[layer_id]                            
        height = self.layer_df.loc[layer_id].height
        radius = math.sqrt(self.layer_df.loc[layer_id].volume/(height/(1000)*math.pi))*1000
        total_neurons = int(round(self.layer_df.loc[layer_id].neuron_count))
        in_ratio = self.layer_df.loc[layer_id].in_ratio
        in_neurons = int(round(total_neurons*(in_ratio/100)))
        ex_neurons = int(round(total_neurons)) - in_neurons
        region = self.layer_df.loc[layer_id].region

        print("Layer: " + str(region))
        print("Radius: " + str(radius))
        print("Height: " + str(height))
        print("Total Neurons: " + str(total_neurons))
        print("IN Neurons: " + str(in_neurons))
        print("EX Neurons: " + str(ex_neurons))

        if(ax is None):
            fig=plt.figure(figsize=(10,10))
            ax = Axes3D(fig, azim=30, elev=30)
                                      
        layer_neuron_df = self.neuron_df.loc[self.neuron_df['layer'] == region]
        for i, row in layer_neuron_df.iterrows():
            neuron_pos = row['pos']
            if(row['type'] == 'EX'):
                ax.scatter(neuron_pos[0], neuron_pos[1], neuron_pos[2], c='r')
            else:
                ax.scatter(neuron_pos[0], neuron_pos[1], neuron_pos[2], c='b')
                                      
        #Create cylindrical layer based on sampled volume and depth
        x = np.linspace(x_center-radius, x_center+radius, resolution)
        z = np.linspace(depth_start, depth_start+height, resolution)
        X, Z = np.meshgrid(x, z)

        Y = np.sqrt(radius**2 - (X - x_center)**2) + y_center # Pythagorean theorem

        ax.plot_surface(X, Y, Z, linewidth=0, color=color, alpha=alpha)
        ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, color=color, alpha=alpha)

        floor = Circle((x_center, y_center), radius, color=color, alpha=alpha)
        ax.add_patch(floor)
        art3d.pathpatch_2d_to_3d(floor, z=depth_start, zdir="z")

        ceiling = Circle((x_center, y_center), radius, color=color, alpha=alpha)
        ax.add_patch(ceiling)
        art3d.pathpatch_2d_to_3d(ceiling, z=depth_start+height, zdir="z")

        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        ax.set_xlabel('x (microns)')
        ax.set_ylabel('y (microns)')
        ax.set_zlabel('Depth (Microns)')
        
    def plotColumn(self, neuron_df=None, plot_neurons=False, resolution=100, alpha=.25, x_center = 0, 
                   y_center = 0, ax=None):
        '''This plots the full column shell and if plot_neurons is True, then it also plots the neurons passed in. If there are
        no neurons passed into the function, but plot_neurons is True, then it will plot all neurons in the column (slow).

        Parameters
        ----------
        neuron_df: DataFrame
            DataFrame containing all the neurons to be plotted.
        plot_neurons: bool
            Boolean that determines whether or not neurons are plotted within the column shell
        resolution: int
            Resolution of the layer shells
        alpha: float
            Transparency of the layer shells
        x_center: float
            x coordinate of the center of the layer shells
        y_center: float
            y coordinate of the center of the layer shells
        ax: Axes3D
            The 3D axes for plotting the column
            
        '''
        radii = []
        heights = []
        neuron_counts = []
        in_counts = []
        ex_counts = []
        layer_names = []
                                      
        for i, row in self.layer_df.iterrows():
            radius = row['radius']
            radii.append(radius)
            height = row['height']
            heights.append(height)    
            region = row['region']
            layer_names.append(region)

        print("Layers: " + str(layer_names))
        print("Radii: " + str(radii))
        print("Heights: " + str(heights))

        if(ax is None):
            fig=plt.figure(figsize=(30,30))
            ax = Axes3D(fig, azim=30, elev=30)

        #Scale axis
        x_scale=2
        y_scale=2
        z_scale=4

        scale=np.diag([x_scale, y_scale, z_scale, 1.0])
        scale=scale*(1.0/scale.max())
        scale[3,3]=1.0

        def short_proj():
            return np.dot(Axes3D.get_proj(ax), scale)

        ax.get_proj=short_proj

        #Plot all neurons
        in_counts = dict((layer,0) for layer in layer_names)
        ex_counts = dict((layer,0) for layer in layer_names)
        total_counts = dict((layer,0) for layer in layer_names)
       
        if(plot_neurons):
            if(neuron_df is None):
                neuron_df = self.neuron_df                                    
            for index, neuron in neuron_df.iterrows():
                neuron_pos = neuron['pos']
                neuron_type = neuron['type']
                neuron_layer = neuron['layer']
                if(neuron_type == 'EX'):
                    color = 'r'
                    ex_counts[neuron_layer] += 1
                    total_counts[neuron_layer] += 1
                else:
                    color = 'b'
                    in_counts[neuron_layer] += 1
                    total_counts[neuron_layer] += 1
                ax.scatter(neuron_pos[0], neuron_pos[1], neuron_pos[2], c=color)
            print("Total Neurons: " + str(total_counts))
            print("IN Neurons: " + str(in_counts))
            print("EX Neurons: " + str(ex_counts))

    #   Plot layer cylinders
        for i in range(len(radii)):
            depth = self.depths[i]
            radius = radii[i]
            height = heights[i]
            #Create cylindrical layer based on sampled volume and depth
            x = np.linspace(x_center-radius, x_center+radius, resolution)
            z = np.linspace(depth, depth+height, resolution)
            X, Z = np.meshgrid(x, z)

            Y = np.sqrt(radius**2 - (X - x_center)**2) + y_center # Pythagorean theorem
            
            color = np.random.rand(3,)
            ax.plot_surface(X, Y, Z, linewidth=0, color=color, alpha=alpha)
            ax.plot_surface(X, (2*y_center-Y), Z, linewidth=0, color=color, alpha=alpha)
            floor = Circle((x_center, y_center), radius, color=color, alpha=alpha)
            ax.add_patch(floor)
            art3d.pathpatch_2d_to_3d(floor, z=depth, zdir="z")

            ceiling = Circle((x_center, y_center), radius, color=color, alpha=alpha)
            ax.add_patch(ceiling)
            art3d.pathpatch_2d_to_3d(ceiling, z=depth+height, zdir="z")

            depth += height
        ax.set_xlim(-300, 300)
        ax.set_ylim(-300, 300)
        ax.set_xlabel('x (microns)')
        ax.set_ylabel('y (microns)')
        ax.set_zlabel('Depth (Microns)')
                                      
    def plotDensityDistribution(self, ax=None):
        '''This plots density vs. depth distribution of the column.

        Parameters
        ----------
        ax: axes
            The 2D axes for plotting the density vs. depth distribution
            
        '''
        
        if(ax is None):
            fig, ax = plt.subplots(ncols=1, nrows=1)
        
        x_ticks = self.depths
        xticks_minor = []
        for i in range(len(self.depths) - 1):
            xticks_minor.append((self.depths[i] + self.depths[i+1])/2)
        xlbls = list(self.layer_df.region.unique())

        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_minor_locator(ticker.FixedLocator(xticks_minor))
        ax.xaxis.set_minor_formatter(ticker.FixedFormatter(xlbls))
   
        ax.set_xticks(self.depths)
        ax.grid(linestyle='--')
        ax.grid(linestyle='--')
        
        plt.plot(self.znew, self.depth_distribution(self.znew), '--')
        plt.xlabel("Depth (microns)")
        plt.ylabel("Density (10^3/mm^3)");
        
    def getNeurons(self):
        '''This returns the neuron DataFrame for the column.

        Returns
        ----------
        self.neuron_df: DataFrame
            DataFrame containing all neurons in the column
            
        '''                              
        return self.neuron_df
    
    def getLayers(self):
        '''This returns the layer DataFrame for the column.

        Returns
        ----------
        self.layer_df: DataFrame
            DataFrame containing all layers in the column
            
        '''                                
        return self.layer_df
    
    def __buildColumn(self, cortex_data_dict):
        '''Private function that fills both the layer Dataframe and neuron DataFrame according to the data dictionary it is 
        passed.
           
           The steps of this function are:
              
              1) Import data_dict for the cortical region
              2) Build density vs depth distributions for each layer using data_dict in order to place neurons later
              3) Fill layer DataFrame using details from data_dict and density/depth distributions
              4) Fill neuron DataFrame using details from data_dict and layer DataFrame

        Parameters
        ----------
        cortex_data_dict: dict
            Contains all biophysical details for creating an average cortical column. Must contain: heights (microns),
            densities (10^3/mm^3), and in_ratios (inhibitory neuron percentage of total).
            
        '''
        #Import data in dictionary format
        attributes = ['heights', 'densities', 'in_ratios']
        new_column_dict = self.__sampleColumnDict(cortex_data_dict, attributes)
        self.num_layers = len(new_column_dict['layers'])
        
        #Generate density vs. depth distributions for the column and for each layer
        print("Building density/depth distributions...")
        self.__generateDensityDistribution(new_column_dict)
        self.__generateLayerDensityDistributions()
        print("... distributions built!")
        
        #Build layer DataFrame
        print("Building layers...")
        layers_list = []
        for layer_id in range(len(new_column_dict['layers'])):
            height = new_column_dict['heights'][layer_id]
            density = new_column_dict['densities'][layer_id]
            volume = (math.pi * (self.radius**2) * height) / (1000.0**3)
            neuron_count = int(((volume) * density)*(10.0**3))
            in_ratio = new_column_dict['in_ratios'][layer_id]
            layer_density_distribution = self.layer_density_distributions[layer_id]
            region = new_column_dict['layers'][layer_id]
            layer = Layer(layer_id, region, height, self.radius, volume, neuron_count, density, in_ratio, 
                          layer_density_distribution)
            layers_list.append(layer)           
        #Layer DataFrame
        self.layer_df = pd.DataFrame(layers_list, columns=Layer._fields)
        print("... layers built!")
        
        #Build neuron DataFrame
        prev_neuron_pos  = None
        all_neurons = []
        all_neuron_count = 0
        
        print("Filling layers with neurons...")
        for layer_id in range(len(new_column_dict['layers'])):
            depth_start = self.depths[layer_id]
            height = new_column_dict['heights'][layer_id]
            total_neurons = self.layer_df[self.layer_df['id'] == layer_id].neuron_count.item()
            in_ratio = new_column_dict['in_ratios'][layer_id]
            in_neurons = int(round(total_neurons*(in_ratio/100)))
            ex_neurons = int(round(total_neurons)) - in_neurons
            total_neurons = in_neurons + ex_neurons
            density_dist = self.layer_density_distributions[layer_id]
            print("...Filling layer " + str(new_column_dict['layers'][layer_id]) + " with " + str(total_neurons) + 
                  " neurons...")
            in_neuron_positions, ex_neuron_positions, neuron_positions = self.__fillNeuronsRejection(in_neurons, ex_neurons, 
                                                                                                     height, depth_start, 
                                                                                                     density_dist,
                                                                                       prev_neuron_positions=prev_neuron_pos) 
            prev_neuron_pos = neuron_positions
            #Fill the neuron dataframe with all the neurons being plotted
            for i in range(in_neuron_positions.shape[0]):
                in_neuron_position = in_neuron_positions[i]
                new_neuron = Neuron(all_neuron_count, in_neuron_position, 'IN', new_column_dict['layers'][layer_id])
                all_neurons.append(new_neuron)
                all_neuron_count += 1
            for i in range(ex_neuron_positions.shape[0]):
                ex_neuron_position = ex_neuron_positions[i]
                new_neuron = Neuron(all_neuron_count, ex_neuron_position, 'EX', new_column_dict['layers'][layer_id])
                all_neurons.append(new_neuron)
                all_neuron_count += 1
        #Neuron DataFrame
        self.neuron_df = pd.DataFrame(all_neurons, columns=Neuron._fields)
        print("... all layers filled!")
    
    def __sampleColumnDict(self, cortex_data_dict, attributes):
        '''Private function for sampling a new column dict given the cortex dict information passed in and the attributes
        needed to define a column. Used as a first step in the buildColumn function. The sampled attributes are the heights 
        (microns), the neurons densities by layer (10^3/mm^3), and the inhibitory to excitatory ratios. 

        Parameters
        ----------
        cortex_data_dict: dict
            Dict containing all information necessary to construct a cortical column.
        attributes: array_like
            List containing all attributes to define a new column.
            
        Returns
        ----------
        sampled_column_dict: dict
            Dict containing biophysical details of a new average column
            
        '''  
        while True:
            sampled_column_dict = {attribute: [] for attribute in attributes}
            valid_column = True
            for attribute in attributes:
                mean = cortex_data_dict['Mean'][attribute]
                cov = np.zeros((len(cortex_data_dict['Mean']['layers']), len(cortex_data_dict['Mean']['layers'])))
                for i, std in enumerate(cortex_data_dict['Std'][attribute]):
                    cov[i][i] = std**2
                sample = np.random.multivariate_normal(mean, cov, 1).T
                sampled_column_dict[attribute] = sample.flatten()

            sampled_column_dict['layers'] = cortex_data_dict['Mean']['layers']

            for attribute in attributes:
                if(sampled_column_dict[attribute][sampled_column_dict[attribute] <= 0].shape[0] > 0):
                    valid_column = False
                if(attribute == 'in_ratios'):
                    if(sampled_column_dict[attribute][sampled_column_dict[attribute] >= 100].shape[0] > 0):
                        valid_column = False
            if valid_column:
                break

        return sampled_column_dict

    def __generateDensityDistribution(self, column_dict):
        '''Private function to generate the density of neurons vs. depth distribution of the column.

        Parameters
        ----------
        column_dict: dict
            Dict containing biophysical details of a average column
            
        '''  
        initial_depth = 0
        distance_discretization = 10
        depths = []
        depth = initial_depth
        depths.append(depth)
        for height in column_dict['heights']:
            depth += height
            depths.append(depth)
        self.depths = depths
                                    
        x = []
        #Scatter/store densities (first and last layer scattered at edges of layer, the rest in middle of layer)
        for i in range(len(column_dict['layers'])):
            if(i > 0 and i < len(column_dict['layers']) - 1):
                x.append((depths[i] + depths[i+1])/2)
            elif(i == 0):
                x.append(depths[i])
            elif(i == len((column_dict['layers'])) - 1):
                x.append(depths[i + 1])

        y = column_dict['densities']
        depth_distribution = interp1d(x, y, kind='cubic')
        znew = np.linspace(x[0], x[-1], num=int(x[-1])*distance_discretization, endpoint=True)
        self.znew = znew
        self.depth_distribution = depth_distribution

    def __generateLayerDensityDistributions(self):
        '''Private function to generate the density of neurons vs. depth distribution for each layer given the 
        depth_distribution.
            
        ''' 
        
        layer_density_distributions = []

        for i in range(self.num_layers):
            densities = self.znew[np.where(np.logical_and(self.znew >= self.depths[i], self.znew < self.depths[i+1]))]
            norm_depth_distribution = [float(depth)/sum(self.depth_distribution(densities)) for depth in self.depth_distribution(densities)]
            xk = np.arange(len(norm_depth_distribution))
            pk = norm_depth_distribution
            layer_density_distribution = stats.rv_discrete(name='custm', values=(xk, pk))
            layer_density_distributions.append(layer_density_distribution)

        self.layer_density_distributions = layer_density_distributions

    def __fillNeuronsRejection(self, in_neurons, ex_neurons, height, depth_start, density_dist, prev_neuron_positions=None):
        '''Private function to fill a column layer with the correct number and type of neuron corresponding to the neuron 
        density distributions for the layer and the min distance between neurons. Rejection sampling based method that uses 
        cKDTrees for speeding up distance search times.

        Parameters
        ----------
        in_neurons: int
            Number of inhibitory neurons that will be placed in the layer
        ex_neurons: int
            Number of excititatory neurons that will be placed in the layer
        height: float
            Height of the layer (microns)
        depth_start: float
            The depth of the base of the layer (microns)
        density_dist: stats.rv_discrete
            The density vs. depth density for the given layer
        prev_neuron_positions: numpy.ndarray
            Array containing all the positions of neurons in the previous layer (to ensure that all newly placed neurons adhere 
            to the min neuron distance even near the boundary of layers). Only None when filling first layer where there are no 
            boundary points about which to consider.
            
        Returns
        ----------
        in_neuron_positions: numpy.ndarray
            Array containing all positions of newly placed inhibitory neurons
        ex_neuron_positions: numpy.ndarray
            Array containing all positions of newly placed excitatory neurons
        curr_neuron_positions: numpy.ndarray
            Array containing all positions of newly placed neurons
            
        '''  
        
        Z_inds = density_dist.rvs(size=in_neurons+ex_neurons)
        Z_ratios = [float(idx) / density_dist.xk.shape[0] for idx in Z_inds]
        Z = np.asarray([Z_ratio * height + depth_start for Z_ratio in Z_ratios])

        #Add first point
        z = Z[0]
        length = np.sqrt(np.random.uniform(0, self.radius**2))
        angle = np.pi * np.random.uniform(0, 2)
        x = length * np.cos(angle)
        y = length * np.sin(angle)
        curr_neuron_positions = np.asarray([x, y, z]).reshape((1,3))

        if(prev_neuron_positions is None):
            all_neuron_positions = None
        else:
            all_neuron_positions = prev_neuron_positions
            all_neuron_positions = np.append(all_neuron_positions, np.asarray([x, y, z]).reshape((1,3)), axis=0)

        #Try to add additional points with min dist requirement between points
        for i, z in enumerate(Z[1:]):
            accepted_sample = False
            while not accepted_sample:
                length = np.sqrt(np.random.uniform(0, self.radius**2))
                angle = np.pi * np.random.uniform(0, 2)
                x = length * np.cos(angle)
                y = length * np.sin(angle)
                if(prev_neuron_positions is None):
                    tree = cKDTree(curr_neuron_positions)
                    idx = tree.query_ball_point((x, y, z), self.min_neuron_dist)
                else:
                    tree = cKDTree(all_neuron_positions)
                    idx = tree.query_ball_point((x, y, z), self.min_neuron_dist)
                if(prev_neuron_positions is None):
                    if(curr_neuron_positions[idx].shape[0] == 0):
                        curr_neuron_positions = np.append(curr_neuron_positions, np.asarray([x, y, z]).reshape((1,3)), axis=0)
                        accepted_sample = True
                else:
                    if(all_neuron_positions[idx].shape[0] == 0):
                        curr_neuron_positions = np.append(curr_neuron_positions, np.asarray([x, y, z]).reshape((1,3)), axis=0)
                        all_neuron_positions = np.append(all_neuron_positions, np.asarray([x, y, z]).reshape((1,3)), axis=0)
                        accepted_sample = True

        curr_neuron_positions = np.asarray(curr_neuron_positions)        
        in_indices = np.random.choice(curr_neuron_positions.shape[0], in_neurons, replace=False)
        in_neuron_positions = curr_neuron_positions[in_indices]
        ex_neuron_positions = np.delete(curr_neuron_positions, in_indices, axis=0)
        
        return in_neuron_positions, ex_neuron_positions, curr_neuron_positions
            
def createRatSomatosensoryCortexDataDict():
    '''Returns dictionary containing all import statistics for recreating an average somatosensory cortical column 
    from a rat. All biophysical details were gathered from these three papers:
                                
                                https://www.ncbi.nlm.nih.gov/pubmed/21949377
                                https://www.ncbi.nlm.nih.gov/pubmed/20534784
                                https://www.ncbi.nlm.nih.gov/pubmed/20534783
    
    Returns
    -------
    cortex_data_dict: dict
        Dictionary containing biophysical details for each layer of a cortical column in a rat somatosensory cortex
    '''
    cortex_data_dict = {'C2':{}, 'D2':{}, 'D3':{}, 'Mean':{}, 'Std':{}}
    layers = ['L1', 'L2', 'L3', 'L4', 'L5A', 'L5B', 'L6A', 'L6B']

    #Layers: L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['C2']['layers'] = layers
    cortex_data_dict['D2']['layers']  = layers
    cortex_data_dict['D3']['layers']  = layers
    cortex_data_dict['Mean']['layers']  = layers
    cortex_data_dict['Std']['layers']  = layers

    #Total densities 10^3/mm^3: L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['C2']['densities'] = [9.6, 99.2, 107.7, 125.0, 54.0, 57.7, 91.9, 37.7]
    cortex_data_dict['D2']['densities'] = [4.9, 81.8, 94.8, 117.1, 50.5, 57.7, 91.7, 43.3]
    cortex_data_dict['D3']['densities'] = [5.2, 78.3, 102.3, 129.6, 58.8, 64.2, 93.3, 44.8]
    cortex_data_dict['Mean']['densities'] = [6.6, 86.4, 101.6, 123.9, 54.4, 59.8, 92.3, 42]
    cortex_data_dict['Std']['densities'] = [2.7, 11.2, 6.5, 6.3, 4.1, 3.8, 0.9, 3.8]

    #Height (μm): L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['C2']['heights'] = [55, 194.2, 186.2, 262.7, 273.9, 273, 295, 171.2]
    cortex_data_dict['D2']['heights'] = [75, 163.9, 298, 294.3, 212.1, 273.2, 284.4, 181.4]
    cortex_data_dict['D3']['heights'] = [95, 154.5, 330.5, 233.1, 216.3, 274.7, 321.9, 203.4]
    cortex_data_dict['Mean']['heights'] = [75, 171, 272, 263, 234, 274, 300, 185]
    cortex_data_dict['Std']['heights'] = [20, 21, 76, 31, 35, 1, 19, 16]
                        
    cortex_data_dict['Mean']['in_ratios'] = [84.3, 17, 9, 8.1, 19.9, 16.2, 8.9, 8.4]
    cortex_data_dict['Std']['in_ratios'] = [9.5, 2.4, 1.1, 0.5, 4.1, 2.1, 0.9, 0.6]
    
#     #NeuN-positive neuron Count: L1, L2, L3, L4, L5A, L5B, L6A, L6B
#     cortex_data_dict['C2']['counts'] = [72, 2619, 2727, 4465, 2011, 2139, 3687, 877]
#     cortex_data_dict['D2']['counts'] = [52, 1897, 4001, 4877, 1517, 2231, 3691, 1113]
#     cortex_data_dict['D3']['counts'] = [65, 1601, 4478, 3999, 1684, 2336, 3980, 1208]
#     cortex_data_dict['Mean']['counts'] = [63, 2039, 3735, 4447, 1737, 2235, 3786, 1066]
#     cortex_data_dict['Std']['counts'] = [10, 524, 905, 439, 251, 99, 168, 170]

#     #Neuron Count Standard Column: L1, L2, L3, L4, L5A, L5B, L6A, L6B
#     cortex_data_dict['C2']['counts_standard'] = [69, 2507, 2610, 4274, 1925, 2047, 3529, 839]
#     cortex_data_dict['D2']['counts_standard'] = [46, 1674, 3530, 4303, 1339, 1969, 3257, 982]
#     cortex_data_dict['D3']['counts_standard'] = [60, 1471, 4115, 3675, 1547, 2147, 3657,1110]
#     cortex_data_dict['Mean']['counts_standard'] = [58, 1884, 3418, 4084, 1604, 2054, 3481, 977]
#     cortex_data_dict['Std']['counts_standard'] = [12, 549, 758, 355, 297, 89, 204, 135]

#     #Volume (mm^3): L1, L2, L3, L4, L5A, L5B, L6A, L6B
#     cortex_data_dict['C2']['volumes'] = [0.007, 0.026, 0.025, 0.036, 0.037, 0.037, 0.040, 0.023]
#     cortex_data_dict['D2']['volumes']  = [0.011, 0.023, 0.042, 0.042, 0.030, 0.039, 0.040, 0.026]
#     cortex_data_dict['D3']['volumes']  = [0.013, 0.020, 0.044, 0.031, 0.029, 0.036, 0.043, 0.027]
#     cortex_data_dict['Mean']['volumes']  = [0.01, 0.023, 0.037, 0.036, 0.032, 0.037, 0.041, 0.025]
#     cortex_data_dict['Std']['volumes']  = [0.003, 0.003, 0.01, 0.005, 0.005, 0.001, 0.001, 0.002]

#     #Excitatory neuron count: L1, L2, L3, L4, L5A, L5B, L6A, L6B, Entire columns: C2, D2, D3)  C2, D2, D3
#     cortex_data_dict['Mean']['ex_counts'] = [10, 1701, 3398, 4089, 1394, 1873, 3449, 976]
#     cortex_data_dict['Std']['ex_counts'] = [7, 484, 807, 425, 247, 53, 165, 155]

#     #IN neuron count: L1, L2, L3, L4, L5A, L5B, L6A, L6B,  C2, D2, D3
#     cortex_data_dict['Mean']['in_counts'] = [53, 338, 338, 358, 343, 362, 337, 90]
#     cortex_data_dict['Std']['in_counts'] = [6, 40, 109, 15, 65, 62, 31, 17]

#     #IN densities 10^3 per mm^3*: L1, L2, L3, L4, L5A, L5B, L6A, L6B)  (C2, D2, D3
#     cortex_data_dict['Mean']['in_densities'] = [5.5, 14.5, 9.1, 10, 10.9, 9.7, 8.2, 3.5]
#     cortex_data_dict['Std']['in_densities'] = [2.2, 0.7, 1.5, 1.1, 3, 1.9, 0.8, 0.4]

    #IN-to-neuron ratios: L1, L2, L3, L4, L5A, L5B, L6A, L6B,  C2, D2, D3
    
    return cortex_data_dict
