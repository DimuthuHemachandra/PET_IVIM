import pandas as pd
import glob
import os

csvs = snakemake.input.avg_csv # Reading avg_csf file paths for all the subjects
algorithms = snakemake.params.algorithms

df_concat_csv_paths= pd.DataFrame()

out_files = []
algos_list = []

def get_df_formatted(df, algorithm, file_name):

    # Pivot the DataFrame
    pivoted_df = df.T

    # Set the first row as column names
    pivoted_df.columns = pivoted_df.iloc[0]
    pivoted_df = pivoted_df[1:]

    # Reset the index
    pivoted_df.reset_index(drop=False, inplace=True)

    pivoted_df.rename(columns={'index': 'subj'}, inplace=True)
    pivoted_df.index.name = None

    pivoted_df.to_csv(f'results/{algorithm}_{file_name}.csv', index=False)

    out_files.append(f'results/{algorithm}_{file_name}.csv')
    algos_list.append(algorithm)



# Loop through each line in the text file
for algorithm in algorithms:

    df_D = pd.DataFrame()
    df_Dstar = pd.DataFrame()
    df_F = pd.DataFrame()

    # Loop through all subjects
    for files in csvs:

        sub_id = files.split(os.path.sep)[-2]

        # Initialize an empty list to hold the DataFrames


        # Read the text file
        df_paths = pd.read_csv(files)

        csv_file_path = df_paths[df_paths['Algorithm']==algorithm]['Mean_metric'].values[0]

        
        # Read the CSV file into a DataFrame and append to the list
        try:
            df = pd.read_csv(csv_file_path)  # Read the CSV file

            df_D['labels'] = df['name'].values
            df_D[sub_id] = df['transformed_D_mean'].values

            df_Dstar['labels'] = df['name'].values
            df_Dstar[sub_id] = df['transformed_Dstar_mean'].values

            df_F['labels'] = df['name'].values
            df_F[sub_id] = df['transformed_F_mean'].values

            print(f"Successfully read: {csv_file_path}")
        except Exception as e:
            print(f"Error reading {csv_file_path}: {e}")

    # Formatting the concatted csvs and saving
    get_df_formatted(df_D, algorithm, 'D')
    get_df_formatted(df_Dstar, algorithm, 'Dstar')
    get_df_formatted(df_F, algorithm, 'F')


df_concat_csv_paths['Algorithm'] = algos_list
df_concat_csv_paths['Concat_file_path'] = out_files

df_concat_csv_paths.to_csv(snakemake.output.concat_csv_files, index=False)