import os
import shutil

# Parent folder path
parent_folder = "pictures_2Repeat"

# Child folder path
acquisition_folder = os.path.join(parent_folder, "acquisition")

# Parameter file path
initial_parameters_file = os.path.join(acquisition_folder, "parameters.txt")

# Destination folder for analysis
analysis_folder = os.path.join(parent_folder, "analysis")

# Iterate through each TIFF image in the acquisition folder to create a parameter file for each image
for filename in os.listdir(acquisition_folder):
    if filename.endswith(".tif"):
        # Get the image name without the extension
        image_name = filename[:-4]
        
        # Check if the parameter file exists for the current image
        image_parameters_file = os.path.join(acquisition_folder, f"{image_name}_parameters.txt")
        if not os.path.isfile(image_parameters_file):
            # Create a new parameter file for the current image
            shutil.copyfile(initial_parameters_file, image_parameters_file)
            
            # Modify the content of the new parameter file
            with open(image_parameters_file, "r+") as file:
                lines = file.readlines()
                modified_lines = []
                for line in lines:
                    if "filename=" in line:
                        # Replace the content after "filename=" with the actual filename
                        modified_line = line.split("filename=")[0] + f"filename={filename}\n"
                        modified_lines.append(modified_line)
                    else:
                        modified_lines.append(line)
                file.seek(0)
                file.writelines(modified_lines)
                file.truncate()

# Iterate through each TIFF image in the acquisition folder to crete a new folder for each image in the analysis folder
for filename in os.listdir(acquisition_folder):
    if filename.endswith(".tif"):
        # Create a new folder for the current image in the analysis folder if it doesn't exist
        image_analysis_folder = os.path.join(analysis_folder, filename)
        # print(image_analysis_folder)
        if not os.path.exists(image_analysis_folder):
            os.makedirs(image_analysis_folder, exist_ok=True)
        kcna1_folder = os.path.join(image_analysis_folder, "kcna1")
        if not os.path.exists(kcna1_folder):
            os.makedirs(kcna1_folder, exist_ok=True)
            print(kcna1_folder)
