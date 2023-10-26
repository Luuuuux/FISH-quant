import os
import pandas as pd
import numpy as np
import bigfish.stack as stack
import bigfish.classification as classification
import bigfish.plot as plot

### Set the path to the image input folder #########################################################
image_path = 'pictures_2Repeat\\acquisition'

### Read in the parameters from the parameter file #################################################
# Get the list of files ending with ".tif" in the folder
original_image_files = [file for file in os.listdir(image_path) if file.endswith('.tif')]
image_files = []
# Add some condition if you want to exclude some files
for item in original_image_files:
        image_files.append(item)
# check if the files are correctly selected (not necessary but useful) 
print(image_files)
print(len(image_files))


# Get initial parameter file
parameters_path = os.path.join(image_path, 'parameters.txt')
print(image_files)
# Iterate over each image file
for image_file in image_files:
    # Create the path to the image file
    new_parameter = image_file[:-4]+'_parameters.txt'
    new_parameter_path = os.path.join(image_path, new_parameter)
    
    # Read the parameters from parameters.txt for the current im

# Read the parameter file
    with open(new_parameter_path, 'r') as file:
        parameters = {}
        for line in file:
            # Skip empty lines and comments
            if line.strip() == '' or line.strip().startswith('#'):
                continue
            
            # Split the line into key and value
            key, value = line.strip().split('=')
            
            # Determine the value's type and convert it if necessary
            if value.strip() == '':
                # Empty value
                converted_value = None
            elif value.strip().isdigit():
                # Integer value
                converted_value = int(value.strip())
            elif value.strip().replace('.', '', 1).isdigit():
                # Float value
                converted_value = float(value.strip())
            elif ',' in value.strip():
                # Tuple value
                converted_value = tuple(map(int, value.strip().split(',')))
            else:
                # String value
                converted_value = value.strip()
            
            # Store the key-value pair in the dictionary
            parameters[key] = converted_value
    
    # Print the dictionary
    print(parameters)
    # load the parameters into the local namespace
    for key, value in parameters.items():
        locals()[key] = value
    #rna_log_filter_sigma_values = [int(x) for x in parameters['rna_log_filter_sigma']]

    # read in the image ############################################################################
    path = os.path.join(path_input, filename)
    image = stack.read_image(path)
    # Check the shape and dtype of the image (useful to debug)
    # usually it should be (z, y, x, c) and uint16
    print("\r shape: {0}".format(image.shape))
    print("\r dtype: {0}".format(image.dtype))
    print(filename)
    # define the output folder #####################################################################
    path_output = path_output_main + "\\" + filename + "\\" + sd_channel_name

    # parse different results files ################################################################
    file_list = os.listdir(path_output)
    original_cell_list = [file for file in file_list if file.endswith(".npz")]
    cell_list = []
    for item in original_cell_list:
        if item.startswith(sd_channel_name):
            cell_list.append(item)
    dataframes = []
    for cell in cell_list:
        # load single cell data
        path = os.path.join(path_output, cell)
        data = stack.read_cell_extracted(path)
        cell_mask = data["cell_mask"]
        nuc_mask = data["nuc_mask"]
        rna_coord = data["rna_coord"]
        foci_coord = data["foci"]
        smfish = data["smfish"]
        
        # compute features
        features, features_names = classification.compute_features(
        cell_mask, nuc_mask, ndim=3, rna_coord=rna_coord,
        smfish=smfish, voxel_size_yx=103,
        foci_coord=foci_coord,
        centrosome_coord=None,
        compute_distance=True,
        compute_intranuclear=True,
        compute_protrusion=True,
        compute_dispersion=True,
        compute_topography=True,
        compute_foci=True,
        compute_area=True,
        return_names=True)
    
        # build dataframe
        features = features.reshape((1, -1))
        df_cell = pd.DataFrame(data=features, columns=features_names)
        dataframes.append(df_cell)
        
    # concatenate dataframes
    df = pd.concat(dataframes)
    
    # save the dataframe in a csv file in the output folder ########################################
    df.reset_index(drop=True, inplace=True)
    path = os.path.join(path_output, sd_channel_name + "_df_features_multicells.csv")
    df.to_csv(path)