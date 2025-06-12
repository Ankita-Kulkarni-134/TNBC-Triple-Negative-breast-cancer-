import os
import subprocess

# Path to your Babel executable
babel_exe = 'babel.exe'

# Input directory containing .sdf files
input_dir = 'Ligands_similarity'  # Specify your directory where .sdf files are located

# Output directory for .smi files
output_dir = 'output_similarity'  # Specify your desired output directory

# Ensure the output directory exists
if not os.path.exists(output_dir):
    os.makedirs(output_dir)


# Function to generate the command for each .sdf file
def run_similarity_search(input_sdf):
    # Extract the base name of the input file (without extension)
    base_name = os.path.splitext(os.path.basename(input_sdf))[0]

    # Define the output .smi file path
    output_smi = os.path.join(output_dir, f'{base_name}.smi')

    # Form the command
    command = [babel_exe, 'natural_unique_drug.fs', output_smi, '-s', input_sdf, '-at83280']

    # Run the command
    subprocess.run(command, check=True)
    print(f"Processed {input_sdf} and saved output as {output_smi}")


# Iterate over all .sdf files in the input directory
for file_name in os.listdir(input_dir):
    if file_name.endswith('.sdf'):
        # Generate the full input file path
        input_sdf = os.path.join(input_dir, file_name)

        # Run the similarity search command for each file
        run_similarity_search(input_sdf)

print("All files processed.")
