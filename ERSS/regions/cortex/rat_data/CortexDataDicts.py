# -*- coding: utf-8 -*-

def createHelmstaedterRatSomatosensoryCortexDataDict():
    '''Returns dictionary containing all import statistics for recreating an average somatosensory cortical column from the Helmsteader a rat. All biophysical details were gathered from these three papers:
                                
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
    
#     #IN densities 10^3 per mm^3*: L1, L2, L3, L4, L5A, L5B, L6A, L6B)  (C2, D2, D3
#     cortex_data_dict['Mean']['in_densities'] = [5.5, 14.5, 9.1, 10, 10.9, 9.7, 8.2, 3.5]
#     cortex_data_dict['Std']['in_densities'] = [2.2, 0.7, 1.5, 1.1, 3, 1.9, 0.8, 0.4]
    
    return cortex_data_dict
