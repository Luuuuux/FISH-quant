import os
import bigfish
import bigfish.stack as stack
import bigfish.multistack as multistack
import bigfish.segmentation as segmentation
import bigfish.plot as plot
print("Big-FISH version: {0}".format(bigfish.__version__))

# Parent folder path
analysis_folder = r"pictures_2Repeat/analysis"

# Iterate through each child folder in the parent folder
for child_folder in os.listdir(analysis_folder):
    child_folder_path = os.path.join(analysis_folder, child_folder)
    child_folder_sd_path = os.path.join(child_folder_path, 'kcna1')
    ### Not necessary to print the child_folder_path
    print( child_folder_path)
    
    # Check if the current item in the parent folder is a directory
    if os.path.isdir(child_folder_path):
        # Check if cell_mask.tif and nuc_mask.tif exist in the child folder
        cell_mask_path = os.path.join(child_folder_path, "cell_mask.tif")
        nuc_mask_path = os.path.join(child_folder_path, "nuc_mask.tif")
        if not (os.path.isfile(cell_mask_path) and os.path.isfile(nuc_mask_path)):
            # Skip the child folder if either of the files is missing
            print(f"Skipping child folder: {child_folder} (missing cell_mask.tif or nuc_mask.tif)")
            continue
        
        print("Processing child folder:", child_folder)
        
        # Read the nuc_mask and cell_mask images
        nuc_mask = stack.read_image(nuc_mask_path)
        cell_mask = stack.read_image(cell_mask_path)
        
        # Perform image processing (traansfromaion from 0,1 to boolean)
        nuc_boolean = nuc_mask.astype(bool)
        cell_boolean = cell_mask.astype(bool)
        # Transfor binary boolean mask to label
        nuc_label = segmentation.label_instances(nuc_boolean)
        cell_label = segmentation.label_instances(cell_boolean)
        # Clean the segmentation (Can adjust the parameters in the function)
        nuc_label = segmentation.clean_segmentation(nuc_label, delimit_instance=True, fill_holes= True, small_object_size = 20)
        cell_label = segmentation.clean_segmentation(cell_label, delimit_instance=True, fill_holes= True, small_object_size = 50)
        # Match the nuc_label and cell_label images
        nuc_label, cell_label = multistack.match_nuc_cell(nuc_label, cell_label, single_nuc=True, cell_alone=False)
        
        # Save the resulting nuc_label and cell_label images
        nuc_label_path = os.path.join(child_folder_path, "nuc_label.tif")
        cell_label_path = os.path.join(child_folder_path, "cell_label.tif")

        nuc_label_sd_path = os.path.join(child_folder_sd_path, "nuc_label.tif")
        cell_label_sd_path = os.path.join(child_folder_sd_path, "cell_label.tif")
        #plot.plot_images([nuc_label, cell_label], titles=["Nucleus", "Cell"])
        stack.save_image(nuc_label, nuc_label_path)
        stack.save_image(cell_label, cell_label_path)
        stack.save_image(nuc_label, nuc_label_sd_path)
        stack.save_image(cell_label, cell_label_sd_path)
        
        print("Processed child folder:", child_folder)

print("Processing completed for all child folders.")
