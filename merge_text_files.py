import os

# Set the directory containing your text files
input_dir = r"E:\TNBC\admet\Rushi\output_files"  # Change this to your directory
output_file = 'merged_output.txt'

# Open the output file in write mode
with open(output_file, 'w', encoding='utf-8') as outfile:
    for filename in os.listdir(input_dir):
        if filename.endswith('.txt'):
            file_path = os.path.join(input_dir, filename)
            with open(file_path, 'r', encoding='utf-8') as infile:
                outfile.write(f"--- {filename} ---\n")  # Optional: write filename as a header
                outfile.write(infile.read())
                outfile.write('\n\n')  # Optional: add space between files

print(f"All .txt files from '{input_dir}' have been merged into '{output_file}'.")
