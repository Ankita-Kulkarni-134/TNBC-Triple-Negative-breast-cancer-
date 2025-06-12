import os

# Define the base filename and split number
f = "C:\Users\adina\OneDrive\Desktop\admet\Rushi\combined_unique.txt"
split_number = 200

# Define the output directory
output_dir = "C:\Users\adina\OneDrive\Desktop\admet\Rushi\output_files"

# Create the output directory if it doesn't exist
if not os.path.exists(output_dir):
    os.makedirs(output_dir)

i = 0
j = 0

# Open the first output file
output_path = os.path.join(output_dir, f + '_' + str(j) + '.txt')
f2 = open(output_path, 'w')

# Process each line in the input file
with open(f + '.txt', 'r') as infile:
    for line in infile:
        f2.write(line)
        i += 1
        if i >= split_number:
            f2.close()
            j += 1
            output_path = os.path.join(output_dir, f + '_' + str(j) + '.txt')
            f2 = open(output_path, 'w')
            i = 0  # Reset the count for the new file

# Close the last file
if f2:
    f2.close()

print(f"Total number of lines: {i + (j * split_number)}")
