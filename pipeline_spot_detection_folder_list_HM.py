import os
import bigfish.stack as stack
import bigfish.detection as detection
import bigfish.plot as plot
import numpy as np

### Set the path to the image input folder #########################################################
image_path = 'picture5hr\\acquisition'
rna_image_list = []
rna_image_mip_list = []   
paths = []
path_outputs = []

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
####
for image_file in image_files:
    # Create the path to the image file (not the parameter file itself)
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
    
    # create the output folder to save the results #################################################
    # you may have different output folders for different channels 

    #
    sd_channel_name ='kcna1'
    #     
    path_output = path_output_main + "\\" + filename + "\\" + sd_channel_name
    path_outputs.append(path_output)

    if not os.path.exists(path_output):
        os.makedirs(path_output)
        print("Folder created successfully.")
    else:
        print("Folder already exists.")
    
    # read in the image and split the channels #####################################################
    path = os.path.join(path_input, filename)
    paths.append(path)

    image = stack.read_image(path)
    # Check the shape and dtype of the image (useful to debug)
    # usually it should be (z, y, x, c) and uint16
    print("\r shape: {0}".format(image.shape))
    print("\r dtype: {0}".format(image.dtype))
    print(filename)
    
    # The image is in the order (z, y, x, c) but we need (c, z, y, x)
    # we need to transpose the image
    image_pre =  np.transpose(image, (3, 0, 1, 2))
    print(image_pre.shape)
    
    # define the channel 
    # REMEMBER: if you want to do spot detection on multiple channels, you need to run this script
    # multiple times
    # Each time you need to change spot_detection_channel to the channel you want to detect
    rna = image_pre[spot_detection_channel -1] # -1 because the index starts from 0
    print("\r shape: {0}".format(rna.shape))
    print("\r dtype: {0}".format(rna.dtype))

    # Preprocessing ################################################################################
    
    image_stretched = stack.rescale(rna)
    print("rna stretched")
    print("\r min: {0}".format(image_stretched.min()))
    print("\r max: {0}".format(image_stretched.max()))
    print("\r shape: {0}".format(image_stretched.shape))
    print("\r dtype: {0}".format(image_stretched.dtype))
    
    # You may wish to plot the projection image to check if the stretching worked
    rna_mip = stack.maximum_projection(rna)
    rna_image_mip_list.append(rna_mip)
    # You may wish to plot the projection image to check if the stretching worked
    # If so, uncomment the following line
    # plot.plot_images(rna_mip, contrast=True, framesize=(5, 5))
    
    # Use rescaled image for spot detection
    rna = image_stretched
    print(rna.shape)
    
    # Do filter on 3D stack
    # The filter type is not limited to log filter
    # Be careful about whether the filter can be applied to 3D image (some of them are only 
    # compatilbe with 2D image)
    # You may change the filter type and paramters here
    rna_background = stack.remove_background_gaussian(rna, sigma=[10, 5, 5])
    rna_logfilter_3D = stack.log_filter(rna, sigma=[2, 1, 1])
    rna_filtered = rna_logfilter_3D
    
    # You may wish to plot the projection image to check if the filter worked
    rna_filtered_mip = stack.maximum_projection(rna_filtered)
    # If so, uncomment the following line
    #plot.plot_images(rna_filtered_mip, contrast=True, framesize=(5, 5))
    
    ## Assign the preprocessed image to rna
    rna = rna_filtered
    rna_image_list.append(rna)
    
# Spot Detection ################################################################################
# Remember to change the parameters in parameters.txt when you are doing analysis in a new set of
# images or you want to change the parameters.
# Convert the object radius in pixel.
spot_radius_px = detection.get_object_radius_pixel(
    voxel_size_nm=sd_voxel_size, 
    object_radius_nm=sd_spot_radius, 
    ndim=3)
print("spot radius (z axis): {:0.3f} pixels".format(spot_radius_px[0]))
print("spot radius (yx plan): {:0.3f} pixels".format(spot_radius_px[-1]))
# Apply LoG filter followed by a Local Maximum algorithm to detect spots in a 2-d or 3-d image.
spots, threshold = detection.detect_spots(
    images=rna_image_list, 
    return_threshold=True, 
    voxel_size= sd_voxel_size,
    spot_radius= sd_spot_radius)
print("detected spots")
print(spots)
print("\r shape: {0}".format(spots[1].shape))
#print("\r dtype: {0}".format(spots.dtype))
print("\r threshold: {0}".format(threshold))
    # Detect dense and bright regions with potential clustered spots and simulate a more realistic
    # number of spots in these regions.
for i in range(len(rna_image_list)):  
    # extract the variables from the list  
    rna_spots = spots[i]
    rna = rna_image_list[i]
    rna_mip = rna_image_mip_list[i]
    path_output = path_outputs[i]
    path = paths[i]
    print(path_output)
    spots_post_decomposition, dense_regions, reference_spot = detection.decompose_dense(
        image=rna, 
        spots=rna_spots, 
        voxel_size=sd_voxel_size, 
        spot_radius=sd_spot_radius, 
        alpha=0.7,  # alpha impacts the number of spots per candidate region
        beta=1,  # beta impacts the number of candidate regions to decompose
        gamma=5)  # gamma the filtering step to denoise the image
    print("detected spots before decomposition")
    print("\r shape: {0}".format(rna_spots.shape))
    print("\r dtype: {0}".format(rna_spots.dtype), "\n")
    print("detected spots after decomposition")
    print("\r shape: {0}".format(spots_post_decomposition.shape))
    print("\r dtype: {0}".format(spots_post_decomposition.dtype))
    
    # Cluster spots and detect relevant aggregated structures.
    spots_post_clustering, clusters = detection.detect_clusters(  ##modifiaction!
        spots=spots_post_decomposition, 
        voxel_size=sd_voxel_size, 
        radius=cluster_radius, 
        nb_min_spots=cluster_nb_min_spots)
    print("detected spots after clustering")
    print("\r shape: {0}".format(spots_post_clustering.shape))
    print("\r dtype: {0}".format(spots_post_clustering.dtype), "\n")
    print("detected clusters")
    print("\r shape: {0}".format(clusters.shape))
    print("\r dtype: {0}".format(clusters.dtype))

    ### output the results #########################################################################
    # save the spots and clusters in png files (use maximum projection image as background)
    path_sd_png = os.path.join(path_output, sd_channel_name + "_sd_result.png")
    plot.plot_detection(rna_mip, 
                        spots=[spots_post_decomposition, clusters[:, :3]], 
                        path_output= path_sd_png,
                        linewidth=[1, 2], 
                        fill=[False, True], 
                        contrast=True, show = False)
    
    # save spots and clusters in csv files
    # There 2 files will be used in the next step (cellular analysis)
    path = os.path.join(path_output, sd_channel_name + "_spots.csv")
    stack.save_data_to_csv(spots_post_clustering, path)
    path = os.path.join(path_output, sd_channel_name +"_clusters.csv")
    stack.save_data_to_csv(clusters, path)
print(threshold)

    
    
    
    