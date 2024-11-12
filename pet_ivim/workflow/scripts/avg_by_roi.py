import numpy as np
import pandas as pd
import nibabel as nib
from nilearn.image import load_img, new_img_like
import os

# Load the label names from the TSV file
label_names_df = pd.read_csv(snakemake.input.label_lookup, sep='\t')

# Load the paths CSV file
paths_df = pd.read_csv(snakemake.input.mean_ivim_csv) 

# Load the mask image (assuming the mask image path remains unchanged)
mask_path = snakemake.input.mask 
mask_img = load_img(mask_path)

# Load the segmented image
segmented_img = load_img(snakemake.input.t1_seg)

out_dir = snakemake.params.out_dir
ivim_work = snakemake.params.ivim_work
ivim_work = os.path.dirname(ivim_work) #dropping the last directory to match with the paths in csvs

# Initialize a dictionary to hold mean results for each algorithm
mean_results = {}

# Function to calculate mean values for an image
def calculate_mean_values(image, mask, segmented_img):
    # Apply the mask to the image
    masked_img = new_img_like(image, np.where(mask.get_fdata() > 0, image.get_fdata(), 0))
    masked_data = masked_img.get_fdata()

    # Get unique labels (excluding background)
    unique_labels = np.unique(segmented_img.get_fdata())
    unique_labels = unique_labels[unique_labels != 0]  # Assuming 0 is background

    # Initialize a dictionary to hold mean values
    mean_values = {}

    # Calculate mean for each label
    for label in unique_labels:
        # Create a mask for the current label
        label_mask = (segmented_img.get_fdata() == label)
        label_values = masked_data[label_mask]

        # Calculate mean for this label
        if label_values.size > 0:
            mean_values[int(label)] = np.mean(label_values)
        else:
            mean_values[int(label)] = np.nan  # Assign NaN if no values
    
    return mean_values

algorithms = np.unique(paths_df['Algorithm'].values)

final_csvs= []
df_paths = pd.DataFrame()

# Loop through each row in the paths DataFrame
for algorithm in algorithms:

    image_paths = paths_df[paths_df['Algorithm']==algorithm]

    mean_dfs = label_names_df
    for image_path in image_paths['Mean_metric'].values:
    
        # Load the specific image
        image = load_img(os.path.join(ivim_work,image_path))

        # Calculate mean values for the current image
        mean_values = calculate_mean_values(image, mask_img, segmented_img)

        # Extract the file name without extension to use as the column name
        base_file_name = os.path.basename(image_path).replace('.nii.gz', '')

        mean_df = pd.DataFrame(list(mean_values.items()), columns=['Label', base_file_name])
        # Merge all mean DataFrames on 'Label'
        mean_dfs = pd.merge(mean_dfs, mean_df, on='Label')


    # Save to CSV
    output_filename = f'{out_dir}_{algorithm}_mean_values.csv'
    mean_dfs.to_csv(output_filename, index=False)
    print(f"Mean values for {algorithm} saved to {output_filename}")
    final_csvs.append(output_filename)

df_paths['Algorithm'] = algorithms
df_paths['Mean_metric'] = final_csvs

df_paths.to_csv(snakemake.output.avg_csv, index=False)
