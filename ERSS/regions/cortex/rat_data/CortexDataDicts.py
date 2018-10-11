# -*- coding: utf-8 -*-

def createHelmstaedterRatSomatosensoryCortexDataDict():
    '''Returns dictionary containing all import statistics for recreating an average somatosensory cortical column from the Helmsteader studies. All biophysical details were gathered from these three papers:

                                https://www.ncbi.nlm.nih.gov/pubmed/21949377
                                https://www.ncbi.nlm.nih.gov/pubmed/20534784
                                https://www.ncbi.nlm.nih.gov/pubmed/20534783

    Returns
    -------
    cortex_data_dict: dict
        Dictionary containing biophysical details for each layer of a cortical column in a rat somatosensory cortex
    '''
    cortex_data_dict = {'Mean':{}, 'Std':{}}
    layers = ['L1', 'L2', 'L3', 'L4', 'L5A', 'L5B', 'L6A', 'L6B']

    #Layers: L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['Mean']['layers']  = layers
    cortex_data_dict['Std']['layers']  = layers

    #Densities (10^3/mm^3): L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['Mean']['densities'] = [6.6, 86.4, 101.6, 123.9, 54.4, 59.8, 92.3, 42]
    cortex_data_dict['Std']['densities'] = [2.7, 11.2, 6.5, 6.3, 4.1, 3.8, 0.9, 3.8]

    #Height (μm): L1, L2, L3, L4, L5A, L5B, L6A, L6B
    cortex_data_dict['Mean']['heights'] = [75, 171, 272, 263, 234, 274, 300, 185]
    cortex_data_dict['Std']['heights'] = [20, 21, 76, 31, 35, 1, 19, 16]

    #IN-to-neuron ratios: L1, L2, L3, L4, L5A, L5B, L6A, L6B,  C2, D2, D3
    cortex_data_dict['Mean']['in_ratios'] = [84.3, 17, 9, 8.1, 19.9, 16.2, 8.9, 8.4]
    cortex_data_dict['Std']['in_ratios'] = [9.5, 2.4, 1.1, 0.5, 4.1, 2.1, 0.9, 0.6]

    cortex_data_dict['majority_in_type'] = [None, None, None, None, None, None, None, None]
    cortex_data_dict['majority_ex_type'] = [None, None, None, None, None, None, None, None]

#     #IN densities 10^3 per mm^3*: L1, L2, L3, L4, L5A, L5B, L6A, L6B)  (C2, D2, D3
#     cortex_data_dict['Mean']['in_densities'] = [5.5, 14.5, 9.1, 10, 10.9, 9.7, 8.2, 3.5]
#     cortex_data_dict['Std']['in_densities'] = [2.2, 0.7, 1.5, 1.1, 3, 1.9, 0.8, 0.4]

    return cortex_data_dict

# -*- coding: utf-8 -*-

def createBlueBrainRatSomatosensoryCortexDataDict():
    '''Returns dictionary containing all import statistics for recreating an average somatosensory cortical column from the Blue Brain studies. All biophysical details were gathered from the Blue Brain portal:

                                https://bbp.epfl.ch/nmc-portal/welcome

    Returns
    -------
    cortex_data_dict: dict
        Dictionary containing biophysical details for each layer of a cortical column in a rat somatosensory cortex
    '''
    cortex_data_dict = {'Mean':{}, 'Std':{}}
    layers = ['L1', 'L2', 'L3', 'L4', 'L5', 'L6']

    #Layers: L1, L2, L3, L4, L5, L6
    cortex_data_dict['Mean']['layers']  = layers
    cortex_data_dict['Std']['layers']  = layers

    #Densities (10^3/mm^3): L1, L2, L3, L4, L5, L6
    cortex_data_dict['Mean']['densities'] = [14.2, 164.6, 83.8, 177.3, 83.9, 131.5]
    cortex_data_dict['Std']['densities'] = [4.2, 18.9, 2.7, 11.9, 8.7, 14.5]

    #Height (μm): L1, L2, L3, L4, L5, L6
    cortex_data_dict['Mean']['heights'] = [165.0, 149.0, 353.0, 190.0, 525.0, 700.0]
    cortex_data_dict['Std']['heights'] = [13.0, 13.0, 14.0, 7.0, 33.0, 48.0]

    #IN-to-neuron ratios: L1, L2, L3, L4, L5, L6 (No std. because only one animal used)
    cortex_data_dict['Mean']['in_ratios'] = [100.0, 15.9, 15.9, 10.3, 17.4, 10.8]
    cortex_data_dict['Std']['in_ratios'] = [0, 0, 0, 0, 0, 0]

    #Majority cell type: L1, L2, L3, L4, L5, L6
    cortex_data_dict['majority_in_type'] = ['L1_HAC', 'L2_LBC', 'L3_LBC', 'L4_LBC', 'L5_MC', 'L6_LBC']
    cortex_data_dict['majority_ex_type'] = [None, 'L2_PC', 'L3_PC', 'L4_PC', 'L5_TTPC1', 'L6_IPC']

    '''
    Explanation of cell type notation: The complete name of a neuron is Lx_PC, Lx_DBC etc. - x
    is a number from 1-6, which denotes the layer. Note that layer 1 contains only inhibitory neurons.

	Excitatory:
	PC 	Pyramidal Cell
	SP 	Star Pyramidal Cell
	SS 	Spiny Stellate Cell
	TTPC1 	Thick-tufted Pyramidal Cell with a late bifurcating apical tuft
	TTPC2 	Thick-tufted Pyramidal Cell with an early bifurcating apical tuft
	UTPC 	Untufted Pyramidal Cell
	STPC 	Slender-tufted Pyramidal Cell
	TPC_L4 	Tufted Pyramidal Cell with apical dendrites terminating in layer 4
	TPC_L1 	Tufted Pyramidal Cell with apical dendrites terminating in layer 1
	IPC 	Pyramidal Cell with inverted apical-like dendrites
	BPC 	Pyramidal Cell with bipolar apical-like dendrites


	Inhibitory:
	DAC 	Descending Axon Cell
	NGC-DA 	Neurogliaform Cell with dense axonal arborization
	NGC-SA 	Neurogliaform Cell with slender axonal arborization
	HAC 	Horizontal Axon Cell
	LAC 	Large Axon Cell
	SAC 	Small Axon Cell
	MC 	Martinotti Cell
	BTC 	Bitufted Cell
	DBC 	Double Bouquet Cell
	BP 	Bipolar Cell
	NGC 	Neurogliaform Cell
	LBC 	Large Basket Cell
	NBC 	Nest Basket Cell
	SBC 	Small Basket Cell
	ChC 	Chandelier Cell
    '''

    return cortex_data_dict
