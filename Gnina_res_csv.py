import os

def extract_remark_data(file_content):
    model_data = []
    model_number = None
    minimizedAffinity = None
    CNNscore = None
    CNNaffinity = None
    
    for line in file_content.split('\n'):
        if line.startswith("MODEL"):
            model_number = int(line.split()[1])
        elif "minimizedAffinity" in line:
            try:
                minimizedAffinity = float(line.split()[-1])
            except ValueError:
                minimizedAffinity = None
        elif "CNNscore" in line:
            try:
                CNNscore = float(line.split()[-1])
            except ValueError:
                CNNscore = None
        elif "CNNaffinity" in line:
            try:
                CNNaffinity = float(line.split()[-1])
            except ValueError:
                CNNaffinity = None
        elif line.startswith("ENDMDL"):
            model_data.append((model_number, minimizedAffinity, CNNscore, CNNaffinity))
            # Reset variables for next model
            model_number = None
            minimizedAffinity = None
            CNNscore = None
            CNNaffinity = None

    return model_data

def process_files(input_folder, output_file):
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(output_file, 'w') as out_f:
        out_f.write("Filename,Model,minimizedAffinity,CNNscore,CNNaffinity\n")
        
        for file_name in os.listdir(input_folder):
            if file_name.endswith('.pdbqt'):
                file_path = os.path.join(input_folder, file_name)
                with open(file_path, 'r') as file:
                    file_content = file.read()
                    model_data = extract_remark_data(file_content)
                    
                    for model in model_data:
                        out_f.write(f"{file_name},{model[0]},{model[1]},{model[2]},{model[3]}\n")

# Update the paths below
input_folder = r"D:\PROJECTS\TNBC\Gnina_pdbqt_out"  # Folder containing GNINA .pdbqt files
output_file = r"D:\PROJECTS\TNBC\Gnina_pdbqt_out\gnina_output_summary.csv"  # Output CSV file

# Run the process
process_files(input_folder, output_file)
