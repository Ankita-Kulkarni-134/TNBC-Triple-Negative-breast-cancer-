import os
import csv
import re

def extract_remark_data(file_content):
    results = []
    model_number = None
    vina_result = None
    
    for line in file_content.split('\n'):
        if line.startswith("MODEL"):
            model_number = int(line.split()[1])
        elif line.startswith("REMARK VINA RESULT:"):
            vina_result = line.strip()
        elif line.startswith("ENDMDL"):
            if vina_result:
                # Extract numerical values from VINA result
                match = re.search(r'REMARK VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', vina_result)
                if match:
                    results.append({
                        'Model': model_number,
                        'VINA_Result_1': match.group(1),
                        'VINA_Result_2': match.group(2),
                        'VINA_Result_3': match.group(3)
                    })
                else:
                    results.append({
                        'Model': model_number,
                        'VINA_Result_1': "No result",
                        'VINA_Result_2': "No result",
                        'VINA_Result_3': "No result"
                    })
                vina_result = None  # Reset for the next model
    
    # Check for data after the last ENDMDL if the file doesn't end with ENDMDL
    if vina_result:
        match = re.search(r'REMARK VINA RESULT:\s+([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)', vina_result)
        if match:
            results.append({
                'Model': model_number,
                'VINA_Result_1': match.group(1),
                'VINA_Result_2': match.group(2),
                'VINA_Result_3': match.group(3)
            })
        else:
            results.append({
                'Model': model_number,
                'VINA_Result_1': "No result",
                'VINA_Result_2': "No result",
                'VINA_Result_3': "No result"
            })
    
    return results

def process_files(input_folder, output_file):
    # Ensure the output directory exists
    output_dir = os.path.dirname(output_file)
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Open the output file for writing
    with open(output_file, 'w', newline='') as out_f:
        writer = csv.writer(out_f)
        writer.writerow(["Filename", "Model", "VINA_Result_1", "VINA_Result_2", "VINA_Result_3"])
        
        # Process each file in the input folder
        for file_name in os.listdir(input_folder):
            if file_name.endswith('.pdbqt'):
                file_path = os.path.join(input_folder, file_name)
                with open(file_path, 'r') as file:
                    file_content = file.read()
                    models_data = extract_remark_data(file_content)
                    
                    for model_data in models_data:
                        writer.writerow([
                            file_name,
                            model_data['Model'],
                            model_data['VINA_Result_1'],
                            model_data['VINA_Result_2'],
                            model_data['VINA_Result_3']
                        ])

# Update these paths as needed
input_folder = r"D:\PROJECTS\TNBC\multi-ligand-docking\output_2"  # Folder containing PDBQT files
output_file = r"D:\PROJECTS\TNBC\multi-ligand-docking\output_2"  # File path for the output CSV

process_files(input_folder, output_file)
