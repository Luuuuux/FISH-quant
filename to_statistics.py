import os
import pandas as pd
import matplotlib.pyplot as plt

### Function: list_to_dataframe ####################################################################
# Description: This function concatenates a list of the data or statistics of a certain group
# (mouse type or condition) into a dataframe
# input: list_to_concat: list of the data or statistics of a certain group (mouse type or condition)
# to concatenate
# output: dataframe: dataframe of the data or statistics of a certain group (mouse type or condition)
def list_to_dataframe(list_to_concat, column_names_list):
    if len(list_to_concat) > 0:
        dataframe = pd.concat(list_to_concat, axis=0)
    else:
        dataframe = pd.DataFrame(columns=column_names_list)
    dataframe.rename(columns={'Image': 'Image'}, inplace=True)
    return dataframe

### Major Function: data_to_statistics #############################################################
# Description: This function reads the summary.csv files of the images in the parent folder and
# calculates statistics of the items in the column names. The data and statistics are saved in a csv
# in the statistics folder. The statistics are also plotted and saved in the plot folder.
# Items: cell_area, nb_rna, nb_transcription_site, nb_foci, avg_spot_ts, avg_spot_foci,
# nb_rna_per_cell_area, nb_rna_in_nuc, nb_rna_in_cell, nb_rna_out_cell, ratio_nb_rna_in_over_out,
# ratio_nb_rna_in_over_all, nb_rna_in_cell_over_all_cell_area, nb_rna_out_cell_over_no_cell_area
# Statistics for each item: mean, median, standard deviation, mean - standard deviation,
# mean + standard deviation
# input:
# analysis_folder: The path to the folder stores the analysis results of the images
# (usually ./analysis)
# statistics_folder: The path to the folder you want to save the data csv and the statistics csv
# (usually ./statistics)
# plot_folder: The path to the folder you want to save the plots of the statistics# (usually ./plot)
# sd_channel_name: The name of the spot detection channel (kcna1, etc.)
# mice1 and mice2: The names of the two mouse types (wt, het, etc.)
# condition1 and condition2: The names of the two conditions (sham, dep, etc.)
# output:
# 12 csv files: The data csv and the statistics csv are saved in the statistics folder 
# (usually ./statistics)
# data csv and statistics csv are for:
# condition1, condition2, mice1_condition1, mice1_condition2, mice2_condition1, mice2_condition2
# 14 plots: The plots of the 14 items are saved in the plot folder (usually ./plot)
def data_to_statistics(analysis_folder,statistics_folder, plot_folder, sd_channel_name, mice1, mice2, condition1, condition2):
    # create lists to store the data and statistics of the 6 groups ################################
    condition1_data = []
    condition2_data = []
    mice1_condition1_data = []
    mice1_condition2_data = []
    mice2_condition1_data = []
    mice2_condition2_data = []

    condition1_stat = []
    condition2_stat = []
    mice1_condition1_stat = []
    mice1_condition2_stat = []
    mice2_condition1_stat = []
    mice2_condition2_stat = []
    
    # list the column names we want
    # you may modify the column names if you want to calculate statistics of other items
    column_names = ['Image', 'Item', 'cell_area', "nb_rna", "nb_transcription_site",
               "nb_foci", "avg_spot_ts", "avg_spot_foci", "nb_rna_per_cell_area",
               "nb_rna_in_nuc", "nb_rna_in_cell", "nb_rna_out_cell",
               "ratio_nb_rna_in_over_out", "ratio_nb_rna_in_over_all",
               "nb_rna_in_cell_over_all_cell_area", "nb_rna_out_cell_over_no_cell_area"]  
         
    # loop through the analysis folder to find the summary.csv files for each image ################
    for root, dirs, files in os.walk(analysis_folder):
        for dir_name in dirs:
            if dir_name.endswith('.tif'):
                # find the path to the summary.csv file
                csv_path = os.path.join(root, dir_name, sd_channel_name, sd_channel_name + '_summary.csv')
                if os.path.isfile(csv_path):
                    # Read the csv file into a dataframe
                    df_data = pd.read_csv(csv_path)

                    # Rename the first column to 'Image' and add the image name to the first column
                    df_data.columns.values[0] = 'Image'
                    image_name = dir_name[:-4]
                    df_data.iloc[:, 0] = image_name

                    # Calculate mean, median, and standard deviation
                    mean_values = df_data.mean()
                    median_values = df_data.median()
                    std_values = df_data.std()
                    mean_minus_sdv_values = mean_values - std_values
                    mean_plus_sdv_values = mean_values + std_values

                    # Create an empty dataframe with the column names
                    df_stat = pd.DataFrame(columns=column_names)

                    row1 = [image_name, 'mean', mean_values['cell_area'], mean_values["nb_rna"],
                            mean_values["nb_transcription_site"], mean_values["nb_foci"],
                            mean_values['avg_spot_ts'], mean_values['avg_spot_foci'],
                            mean_values['nb_rna_per_cell_area'], mean_values['nb_rna_in_nuc'],
                            mean_values["nb_rna_in_cell"], mean_values["nb_rna_out_cell"],
                            mean_values["ratio_nb_rna_in_over_out"], mean_values["ratio_nb_rna_in_over_all"],
                            mean_values['nb_rna_in_cell_over_all_cell_area'], mean_values['nb_rna_out_cell_over_no_cell_area']
                        ]
                    row2 = [
                            image_name, 'median', median_values['cell_area'], median_values['nb_rna'],
                            median_values['nb_transcription_site'], median_values['nb_foci'],
                            median_values['avg_spot_ts'], median_values['avg_spot_foci'],
                            median_values['nb_rna_per_cell_area'], median_values['nb_rna_in_nuc'],
                            median_values["nb_rna_in_cell"], median_values["nb_rna_out_cell"],
                            median_values["ratio_nb_rna_in_over_out"], median_values["ratio_nb_rna_in_over_all"],
                            median_values['nb_rna_in_cell_over_all_cell_area'], median_values['nb_rna_out_cell_over_no_cell_area']
                        ]   
                    row3 = [
                        image_name, 'std', std_values['cell_area'], std_values['nb_rna'],
                        std_values['nb_transcription_site'], std_values['nb_foci'],
                        std_values['avg_spot_ts'], std_values['avg_spot_foci'],
                        std_values['nb_rna_per_cell_area'], std_values['nb_rna_in_nuc'],
                        std_values["nb_rna_in_cell"], std_values["nb_rna_out_cell"],
                        std_values["ratio_nb_rna_in_over_out"], std_values["ratio_nb_rna_in_over_all"],
                        std_values['nb_rna_in_cell_over_all_cell_area'], std_values['nb_rna_out_cell_over_no_cell_area']
                    ]
                    row4 = [
                        image_name, 'mean-srd', mean_minus_sdv_values['cell_area'],
                        mean_minus_sdv_values['nb_rna'], mean_minus_sdv_values['nb_transcription_site'],
                        mean_minus_sdv_values['nb_foci'], mean_minus_sdv_values['avg_spot_ts'],
                        mean_minus_sdv_values['avg_spot_foci'], mean_minus_sdv_values['nb_rna_per_cell_area'],
                        mean_minus_sdv_values['nb_rna_in_nuc'], mean_minus_sdv_values["nb_rna_in_cell"],
                        mean_minus_sdv_values["nb_rna_out_cell"], mean_minus_sdv_values["ratio_nb_rna_in_over_out"],
                        mean_minus_sdv_values["ratio_nb_rna_in_over_all"],
                        mean_minus_sdv_values['nb_rna_in_cell_over_all_cell_area'],
                        mean_minus_sdv_values['nb_rna_out_cell_over_no_cell_area']
                    ]
                    row5 = [
                        image_name, 'mean+std', mean_plus_sdv_values['cell_area'], mean_plus_sdv_values['nb_rna'],
                        mean_plus_sdv_values['nb_transcription_site'], mean_plus_sdv_values['nb_foci'],
                        mean_plus_sdv_values['avg_spot_ts'], mean_plus_sdv_values['avg_spot_foci'],
                        mean_plus_sdv_values['nb_rna_per_cell_area'], mean_plus_sdv_values['nb_rna_in_nuc'],
                        mean_plus_sdv_values["nb_rna_in_cell"], mean_plus_sdv_values["nb_rna_out_cell"],
                        mean_plus_sdv_values["ratio_nb_rna_in_over_out"], mean_plus_sdv_values["ratio_nb_rna_in_over_all"],
                        mean_plus_sdv_values['nb_rna_in_cell_over_all_cell_area'], mean_plus_sdv_values['nb_rna_out_cell_over_no_cell_area']
                    ]

                    # Add the rows to the dataframe
                    df_stat.loc[len(df_stat)] = row1
                    df_stat.loc[len(df_stat)] = row2
                    df_stat.loc[len(df_stat)] = row3
                    df_stat.loc[len(df_stat)] = row4
                    df_stat.loc[len(df_stat)] = row5

                    # Add the dataframe to the list based on the mice types and condition
                    if condition2 in dir_name:
                        condition2_data.append(df_data)
                        condition2_stat.append(df_stat)
                    if condition1 in dir_name:
                        condition1_data.append(df_data)
                        condition1_stat.append(df_stat)
                    if mice1 in dir_name and condition1 in dir_name:
                        mice1_condition1_data.append(df_data)
                        mice1_condition1_stat.append(df_stat)
                    if mice1 in dir_name and condition2 in dir_name:
                        mice1_condition2_data.append(df_data)
                        mice1_condition2_stat.append(df_stat)
                    if mice2 in dir_name and condition1 in dir_name:
                        mice2_condition1_data.append(df_data)
                        mice2_condition1_stat.append(df_stat)
                    if mice2 in dir_name and condition2 in dir_name:
                        mice2_condition2_data.append(df_data)
                        mice2_condition2_stat.append(df_stat)
                    print('Done with the image: ' + image_name)
                    print('-' * 100)
    
    ### Create the dataframes for the different conditions and mice types ##########################
    condition1_df = list_to_dataframe(condition1_data, column_names_list = column_names)
    condition2_df = list_to_dataframe(condition2_data, column_names_list = column_names)
    mice1_condition1_df = list_to_dataframe(mice1_condition1_data, column_names_list = column_names)
    mice1_condition2_df = list_to_dataframe(mice1_condition2_data, column_names_list = column_names)
    mice2_condition1_df = list_to_dataframe(mice2_condition1_data, column_names_list = column_names)
    mice2_condition2_df = list_to_dataframe(mice2_condition2_data, column_names_list = column_names)

    condition1_df_stat = list_to_dataframe(condition1_stat, column_names_list = column_names)
    condition2_df_stat = list_to_dataframe(condition2_stat, column_names_list = column_names)
    mice1_condition1_df_stat = list_to_dataframe(mice1_condition1_stat, column_names_list = column_names)
    mice1_condition2_df_stat = list_to_dataframe(mice1_condition2_stat, column_names_list = column_names)
    mice2_condition1_df_stat = list_to_dataframe(mice2_condition1_stat, column_names_list = column_names)
    mice2_condition2_df_stat = list_to_dataframe(mice2_condition2_stat, column_names_list = column_names)

    ### Save the dataframes to csv files in statistics folder#######################################
    condition1_df.to_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition1 + '_data.csv'), index=False)
    condition2_df.to_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition2 + '_data.csv'), index=False)
    mice1_condition1_df.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition1 + '_data.csv'), index=False)
    mice1_condition2_df.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition2 + '_data.csv'), index=False)
    mice2_condition1_df.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition1 + '_data.csv'), index=False)
    mice2_condition2_df.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition2 + '_data.csv'), index=False)

    condition1_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition1 + '_stat.csv'), index=False)
    condition2_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition2 + '_stat.csv'), index=False)
    mice1_condition1_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition1 + '_stat.csv'), index=False)
    mice1_condition2_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition2 + '_stat.csv'), index=False)
    mice2_condition1_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition1 + '_stat.csv'), index=False)
    mice2_condition2_df_stat.to_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition2 + '_stat.csv'), index=False)
    
    print('Done with the csv files! You may check the statistics folder to see the results later.')
    print('-' * 100)
    
    ### Plot: List of items to plot folder #########################################################
    items = ['cell_area', "nb_rna", "nb_transcription_site",
            "nb_foci", "avg_spot_ts", "avg_spot_foci", "nb_rna_per_cell_area",
            "nb_rna_in_nuc", "nb_rna_in_cell", "nb_rna_out_cell",
            "ratio_nb_rna_in_over_out", "ratio_nb_rna_in_over_all",
            "nb_rna_in_cell_over_all_cell_area", "nb_rna_out_cell_over_no_cell_area"]

    # Iterate over items and create box plots
    for i, item in enumerate(items):
        # Create a new figure for each item
        fig, ax = plt.subplots()
        
        # Extract the mean values for the different groups
        mean_values_condition1 = condition1_df_stat.loc[condition1_df_stat['Item'] == 'mean', item].values
        mean_values_condition2 = condition2_df_stat.loc[condition2_df_stat['Item'] == 'mean', item].values
        mean_values_mice1_condition1 = mice1_condition1_df_stat.loc[mice1_condition1_df_stat['Item'] == 'mean', item].values
        mean_values_mice1_condition2 = mice1_condition2_df_stat.loc[mice1_condition2_df_stat['Item'] == 'mean', item].values
        mean_values_mice2_condition1 = mice2_condition1_df_stat.loc[mice2_condition1_df_stat['Item'] == 'mean', item].values
        mean_values_mice2_condition2 = mice2_condition2_df_stat.loc[mice2_condition2_df_stat['Item'] == 'mean', item].values
        
        # Create a box plot for the current statistic
        ax.boxplot([mean_values_condition1, mean_values_condition2, mean_values_mice1_condition1, mean_values_mice1_condition2,
                    mean_values_mice2_condition1, mean_values_mice2_condition2])
        ax.set_xticklabels([condition1, condition2, mice1 + '_' + condition1, mice1 + '_' + condition2, mice2 + '_' + condition1, mice2 + '_' + condition2], rotation = 45)
        ax.set_title(item)
        ax.set_xlabel('Group')
        ax.set_ylabel('Values')
    
        # Adjust the layout and spacing
        plt.tight_layout()
        
        # Save the plot to plot folder
        plt.savefig(os.path.join(plot_folder, f'{sd_channel_name}_{item}.png'))
        
        # Show the plot if you want (just uncomment the line below)
        # plt.show()

    print('Done with the plots! You may check the plot folder to see the results later.')
    print('-' * 100)

### Run the script! ################################################################################
if __name__ == "__main__":
    print('This script is being run directly.')
    analysis_folder = r'E:\My Drive\pictures_2Repeat\analysis'
    statistics_folder = r'E:\My Drive\pictures_2Repeat\statistics'
    plot_folder = r'E:\My Drive\pictures_2Repeat\plot'
    sd_channel_name = "kcna1"

    # Call the function to process folders, generate CSV files and draw plots
    data_to_statistics(analysis_folder, statistics_folder, plot_folder, sd_channel_name,
                       mice1 = 'wt',mice2 = 'het',
                       condition1 = 'sham', condition2 ='dep')