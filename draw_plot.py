import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import os

### Function: mean_and_stdto #######################################################################
# Description: This function calculates the mean and standard deviation of a given column in a dataframe
# Input: 
# dataframe: the dataframe to calculate the mean and standard deviation from
# stat: the column name to calculate the mean and standard deviation from
# Output:
# mean: the mean of the column
# std: the standard deviation of the column
# !!! Note: if the dataframe is empty, the mean and standard deviation will be 0
def mean_and_std(dataframe, item):
    if dataframe.empty:
        return 0, 0
    mean = dataframe[item].mean()
    std = dataframe[item].std()
    return mean, std

### Function: error_bar ############################################################################
# Description: This function adds mean and mean ± standard deviation as error bars to a given axis
# Input:
# order: the order of the groups in the plot
# group: the group to add the error bars to
# ax: the axis to add the error bars to
# group_means: the mean of the group
# group_std: the standard deviation of the group
# !!! Note: the empty groups will only have a horizontal line at the 0 
# (since the mean and standard deviation are 0)
def error_bar(order, group, ax, group_means, group_std):
    group_pos = order.index(group)
    ax.errorbar([group_pos], [group_means], yerr=[group_std], color="black", fmt="-")
    ax.hlines(y=group_means, xmin=group_pos - 0.1, xmax=group_pos + 0.1, color="black", linewidth=1)
    ax.hlines(y=[group_means - group_std, group_means + group_std], xmin=group_pos - 0.1, xmax=group_pos + 0.1, color="black", linewidth=1)

### Major Function: scatter ########################################################################
# Description: This function creates a scatter plot with error bars for the given statistics
# Input:
# statistics_folder: the folder containing the data csv files we genetaed in the previous step
# plot_folder: the folder to save the plots to
# sd_channel_name: the name of the spot detection channel
# mice1, mice 2 : the types of the mice
# condition1, condition2: the conditions of the mice
# Remarks: The plot is generated with the package seaborn (scatter plot)
# You may check the image each scatter spots belongs to by the spot color and the legend of the plot
def scatter(statistics_folder, plot_folder, sd_channel_name,mice1, mice2, condition1, condition2):
    ### read the data csv files to create dataframes ###############################################
    condition1_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition1 + '_data.csv'))
    condition2_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name + '_' + condition2 + '_data.csv'))
    mice1_condition1_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition1 + '_data.csv'))
    mice1_condition2_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice1 + '_' + condition2 + '_data.csv'))
    mice2_condition1_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition1 + '_data.csv'))
    mice2_condition2_df = pd.read_csv(os.path.join(statistics_folder, sd_channel_name +'_' + mice2 + '_' + condition2 + '_data.csv'))
    
    # get the column names of the dataframes
    selected_columns = condition1_df.iloc[:, 2:]
    column_names = selected_columns.columns.tolist()

    # print to check the column names
    print(column_names)
     
    ### combine the dataframes to create a new csv with data from all the 4 groups #################
    combined_df = pd.concat([mice1_condition1_df, mice1_condition2_df, mice2_condition1_df, mice2_condition2_df])

    # Add a new column "Group" based on the file names
    combined_df["Group"] = combined_df["Image"].apply(lambda x: (mice1 + "_" + condition1) if (mice1 in x) and (condition1 in x) else
                                             (mice1 + "_" + condition2) if (mice1 in x) and (condition2 in x) else
                                             (mice2 + "_" + condition1) if (mice2 in x) and (condition1 in x) else
                                             (mice2 + "_" + condition2) if (mice2 in x) and (condition2 in x) else "")
    
    # output the combined dataframe to a csv file in the statistics folder
    combined_df.to_csv(statistics_folder + "\\" + sd_channel_name + '_combined_plot.csv', index=False)
    print('The combined dataframe is saved to ' + statistics_folder + "\\" + sd_channel_name + '_combined_plot.csv')
    print('-' * 100)
    
    ### draw the scatter plot ######################################################################
    for item in column_names:
        # calculate the mean and standard deviation of each group
        mice1_condition1_means, mice1_condition1_std = mean_and_std(mice1_condition1_df, item)
        mice1_condition2_means, mice1_condition2_std = mean_and_std(mice1_condition2_df, item)
        mice2_condition1_means, mice2_condition1_std = mean_and_std(mice2_condition1_df, item)
        mice2_condition2_means, mice2_condition2_std = mean_and_std(mice2_condition2_df, item)

        # draw the frame of scatter plot with the order in the list 'order'
        sns.set(style="ticks")
        order = [(mice1 + "_" + condition1), (mice1 + "_" + condition2), (mice2 + "_" + condition1), (mice2 + "_" + condition2)]
        scatter_plot = sns.catplot(x = "Group", y = item, hue="Image", data = combined_df.reset_index(), kind = "swarm", order = order)

        #Add mean and mean ± standard deviation as error bars to the scatter plot
        for ax in scatter_plot.axes.flatten():
            error_bar(order, (mice1 + "_" + condition1), ax, mice1_condition1_means, mice1_condition1_std)
            error_bar(order, (mice1 + "_" + condition2), ax, mice1_condition2_means, mice1_condition2_std)
            error_bar(order, (mice2 + "_" + condition1), ax, mice2_condition1_means, mice2_condition1_std)
            error_bar(order, (mice2 + "_" + condition2), ax, mice2_condition2_means, mice2_condition2_std)

        plt.subplots_adjust(top=0.9,bottom = 0.1)
        plt.suptitle(sd_channel_name +"_" + item)
        plt.xlabel("Group")
        plt.ylabel(item)

        # Save the plot in SVG format
        plt.savefig(os.path.join(plot_folder, sd_channel_name + "_" + item + '.svg'), format='svg')
        
        # show the plot if you want (just uncomment the line below)
        # plt.show() 
        print('The plot for ' + sd_channel_name + "_" + item + ' is saved to ' + plot_folder)
        print('-' * 100)

if __name__ == "__main__":
    print('This script is being run directly.')
    statistics_folder = r'E:\My Drive\pictures_2Repeat\statistics'
    plot_folder = r'E:\My Drive\pictures_2Repeat\plot'
    sd_channel_name = "kcna1"
    scatter(statistics_folder, plot_folder, sd_channel_name, 
            mice1 = 'wt', mice2 = 'het',
            condition1 = 'sham',condition2 = 'dep')