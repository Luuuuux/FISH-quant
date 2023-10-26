import os
import numpy as np
import bigfish
import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.plot as plot
import pandas as pd

### Set the path to the image input folder #########################################################
image_path = 'pictures_2Repeat\\acquisition'

### !!! Close the spots.cvs and clusters.csv BEFORE running this script !!! ########################
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
# Iterate over each image file
for image_file in image_files:
    # Create the path to the image file
    new_parameter = image_file[:-4]+'_parameters.txt'
    new_parameter_path = os.path.join(image_path, new_parameter)
    
    # Read the parameters from parameters.txt for the current image
    with open(new_parameter_path, 'r') as file:
        parameters = {}
        # Iterate over each line in the file
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

    # grab the parameters from the dictionary
    for key, value in parameters.items():
        locals()[key] = value
    # special parameters : some functions may need multiple parameters like (,,) or [,,]
    # we grab this kind of parameters by tuple (,,).
    # if you have some parameters used in [,,], you need to convert it to list.
    rna_log_filter_sigma_values = [int(x) for x in parameters['rna_log_filter_sigma']]

    # read in the image and split the channels #####################################################
    path = os.path.join(path_input, filename)
    image = stack.read_image(path)
    # Check the shape and dtype of the image (useful to debug)
    # usually it should be (z, y, x, c) and uint16
    print("\r shape: {0}".format(image.shape))
    print("\r dtype: {0}".format(image.dtype))
    print(filename)
    
    # The image is in the order (z, y, x, c) but we need (c, z, y, x)
    image_pre =  np.transpose(image, (3, 0, 1, 2))
    #print(image_pre.shape)
    #channel1 = image_pre[0]
    #channel2 = image_pre[1]
    #channel3 = image_pre[2]
    #image_list = [channel1,channel2,channel3]
    #print("\r shape: {0}".format(image.shape))
    
    rna = image_pre[rna_channel - 1] # rna channel
    nuc = image_pre[nuc_channel - 1]
    cell = image_pre[cell_channel - 1]

    # preprocess the image and project each channel ################################################
    # preprocess rna channel to make it easy to see the spots in the projection image
    # may change the filetring parameters and methods
    rna_filtered = stack.log_filter(rna, sigma=[2.25,1.5,1.5])
    rna_mip = stack.maximum_projection(rna_filtered)
    rna_mip_stretched = stack.rescale(rna_mip, channel_to_stretch=0)
    # You may want to check the projection image
    # if so, uncomment the following line
    #plot.plot_images(rna_mip_stretched, contrast=True, framesize=(5, 5))
    
    # project the other channels
    nuc_mip = stack.maximum_projection(nuc)
    cell_mip = stack.maximum_projection(cell)

    # assign rna_mip_stretched to image_contrasted, which will be used as the background image in
    # in our results
    # you can change the background image to other channels if needed
    image_contrasted = rna_mip_stretched
    
    # segmented cells
    ## need to adjust when change data
    path_output = path_output_main + "\\" + filename + "\\" + sd_channel_name
    path_output_pic = os.path.join(path_output, sd_channel_name + "_" + image_file[:-4] + "_cellular_results")

    # Input : 1. spot detection results (spots.csv clusters.csv) ###################################
    ########  2. segmented cell and nuclei images (cell_label.tif nuc_label.tif) ###################
    # from the output folder (output of segmentation and spot detetcion)
     
    # segmented cells
    path = os.path.join(path_output, "cell_label.tif")
    cell_label = stack.read_image(path)
    print("segmented cells")
    print("\r shape: {0}".format(cell_label.shape))
    print("\r dtype: {0}".format(cell_label.dtype), "\n")
    
    # segmented nuclei
    path = os.path.join(path_output, "nuc_label.tif")
    nuc_label = stack.read_image(path)
    print("segmented nuclei")
    print("\r shape: {0}".format(nuc_label.shape))
    print("\r dtype: {0}".format(nuc_label.dtype), "\n")
    
    # detected spots
    path = os.path.join(path_output, sd_channel_name +"_spots.csv")
    spots = stack.read_array_from_csv(path, dtype=np.int64)
    print("detected spots")
    print("\r shape: {0}".format(spots.shape))
    print("\r dtype: {0}".format(spots.dtype), "\n")
    
    # detected foci
    path = os.path.join(path_output, sd_channel_name +"_clusters.csv")
    clusters = stack.read_array_from_csv(path, dtype=np.int64)
    print("detected clusters")
    print("\r shape: {0}".format(clusters.shape))
    print("\r dtype: {0}".format(clusters.dtype))
    
    # transcripion sites and foci assignment #######################################################

    # Sometimes, if the spot detection result is too sparse, there will be no (or only one) foci
    # detected. If we do not do some special treatment, some functions will raise an error.

    # No cluster detected
    # We need to assign an empty array to the clusters.
    if len(clusters) == 0:
        clusters =  np.empty((0, 5),dtype= np.int64)

    # Only one cluster detected
    # In this case, the shape of the clusters will be (5,). We need to convert it to (1, 5).
    # go back to check the array shape if you have some errors!!!
    if clusters.shape == (5,):
        clusters = np.reshape(clusters, (1, 5))

    # assign transcription sites and foci ##########################################################
    # ts (transciption sites) = clusters inside nuclei
    # foci = clusters outside nuclei but inside cells
    spots_no_ts, foci, ts = multistack.remove_transcription_site(spots, clusters, nuc_label, ndim=3)
    print("detected spots (without transcription sites)")
    print("\r shape: {0}".format(spots_no_ts.shape))
    print("\r dtype: {0}".format(spots_no_ts.dtype))
    
    spots_in, spots_out = multistack.identify_objects_in_region(nuc_label, spots, ndim=3)
    print("detected spots (inside nuclei)")
    print("\r shape: {0}".format(spots_in.shape))
    print("\r dtype: {0}".format(spots_in.dtype), "\n")
    print("detected spots (outside nuclei)")
    print("\r shape: {0}".format(spots_out.shape))
    print("\r dtype: {0}".format(spots_out.dtype))
    
    # extract cell features ########################################################################
    fov_results = multistack.extract_cell(
        cell_label=cell_label, 
        ndim=3, 
        nuc_label=nuc_label, 
        rna_coord=spots_no_ts, 
        others_coord={"foci": foci, "transcription_site": ts},
        image=image_contrasted,
        others_image={"dapi": nuc_mip, "smfish": rna_mip})
    print("number of cells identified: {0}".format(len(fov_results)))
    
    # create lists to save the average rna intensity in each ts and foci
    ts_avg_rna = []
    foci_avg_rna = []
    
    # create cell_results dataframe to store all the cell features
    for i, cell_results in enumerate(fov_results):
        print("cell {0}".format(i))
        # store cell results
        cell_mask = cell_results["cell_mask"]
        cell_coord = cell_results["cell_coord"]
        # store nucleus results
        nuc_mask = cell_results["nuc_mask"]
        nuc_coord = cell_results["nuc_coord"]
        # store rna results
        rna_coord = cell_results["rna_coord"]
        # store ts and foci results
        ts_coord = cell_results["transcription_site"]
        foci_coord = cell_results["foci"]
        # store image results
        image_contrasted = cell_results["image"]

        # save ts and foci coordinates of each ts and foci in lists
        ts = ts_coord.tolist()
        foci = foci_coord.tolist()
        print(ts)
        print(foci)
        
        # Calculate the mean of number of rna of each ts and foci
        # If no ts or foci are detected (len(ts) == 0 or len(foci) == 0),
        # the mean intensity will be 0
        # Fisrst check the number of rna, ts and foci detected
        print("\r number of rna {0}".format(len(rna_coord)))
        print("\r number of foci {0}".format(len(foci_coord)))
        print("\r number of transcription sites {0}".format(len(ts_coord)))

        if len(ts) > 0:
            # extract the fourth element of each ts sublist
            fourth_elements_ts = [sublist[3] for sublist in ts]
            print(fourth_elements_ts)
            mean_ts = sum(fourth_elements_ts)/len(fourth_elements_ts)
            print(mean_ts)
            ts_avg_rna.append(mean_ts)
        else:
            ts_avg_rna.append(0)

        if len(foci) > 0:
            # extract the fourth element of each foci sublist
            fourth_elements_foci = [sublist[3] for sublist in foci]
            print(fourth_elements_foci)
            mean_foci = sum(fourth_elements_foci)/len(fourth_elements_foci)
            print(mean_foci)
            foci_avg_rna.append(mean_foci)
        else:
            foci_avg_rna.append(0)
                
    # plot cell and save the results of each cell in a png image ###################################
        plot.plot_cell(
           ndim=3, boundary_size = 2, cell_coord=cell_coord, nuc_coord=nuc_coord, 
           rna_coord=rna_coord, foci_coord=foci_coord, other_coord=ts_coord, 
           image=image_contrasted, cell_mask=cell_mask, nuc_mask=nuc_mask, 
           title="Cell {0}".format(i),show = False, path_output = path_output_pic + "_cell" + str(i) + '.png')
    

    # save the cell results in a dataframe called df ###############################################
    # Meaning of the columns:

    # Columns for Single Cells:
    # cell_id: cell id
    # cell_area: cell area (in pixels)
    # nuc_area: nucleus area (in pixels)
    # nb_rna: number of rna in the cell
    # nb_rna_in_nuc: number of rna in the nucleus
    # nb_rna_out_nuc: number of rna in the cytoplasm
    # nb_foci: number of foci in the cell
    # nb_transcription_site: number of transcription sites in the cell
    # avg_spot_foci: average number of rna in each foci
    # avg_spot_ts: average number of rna in each transcription site
    # nb_rna_per_cell_area: number of rna per cell area (in pixels)

    # Columns for a certain image (whole image) (not for single cells):
    # cell_out_area: all the area outside the cells (in pixels)
    # nb_rna_in_cell: number of rna in ALL the cells
    # nb_rna_out_cell: number of rna outside ALL the cells
    # ratio_nb_rna_in_over_out: ratio of number of rna in ALL the cells over the number of rna out of 
    # ALL the cells
    # ratio_nb_rna_in_over_all: ratio of number of rna in ALL the cells over the number of rna in the 
    # whole image
    # nb_rna_in_cell_over_all_cell_area: number of rna in ALL the cells over the cell area (in pixels)
    # nb_rna_out_cell_over_cell_out_area: number of rna out of ALL the cells over the area
    # out of ALL the cell area (in pixels)

    df = multistack.summarize_extraction_results(fov_results, ndim=3)
    print(ts_avg_rna)
    print(foci_avg_rna)
    df["avg_spot_foci"] = foci_avg_rna
    df['avg_spot_ts'] = ts_avg_rna
    df['nb_rna_per_cell_area'] = df['nb_rna'] / df['cell_area']
    nb_rna_all = spots.shape[0]
    nb_rna_in_cell = df['nb_rna'].sum()
    nb_rna_out_cell = nb_rna_all - nb_rna_in_cell
    area_all = cell_label.shape[0] * cell_label.shape[1] 
    cell_area_all = df['cell_area'].sum()
    no_cell_area = area_all - cell_area_all
    df['cell_out_area'] = area_all - cell_area_all
    df['nb_rna_in_cell'] = nb_rna_in_cell
    df['nb_rna_out_cell'] = nb_rna_out_cell
    df["ratio_nb_rna_in_over_out"] = nb_rna_in_cell/nb_rna_out_cell
    df["ratio_nb_rna_in_over_all"] = nb_rna_in_cell/nb_rna_all
    df['nb_rna_in_cell_over_all_cell_area'] = nb_rna_in_cell/cell_area_all
    df['nb_rna_out_cell_over_no_cell_area'] = nb_rna_out_cell/no_cell_area
    print("shape: {0}".format(df.shape))
    print(df.head())
    
    # save the df in a csv file with the name of the channel of spot detection and the word summary
    path = os.path.join(path_output, sd_channel_name + "_summary.csv")
    df.to_csv(path)

    # save the df in a npz file with the name of the channel of spot detection
    for i, cell_results in enumerate(fov_results):
        # save results of each cell in a npz file
        path = os.path.join(path_output, sd_channel_name +"_results_cell_{0}.npz".format(i))
        stack.save_cell_extracted(cell_results, path)
    
    
    
    
    
    