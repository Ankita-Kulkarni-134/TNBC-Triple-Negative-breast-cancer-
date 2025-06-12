
import os
import requests

# Define the save folder
save_folder = r"D:\PROJECT MSC BIOINFO\TNBC\impactdrugs"
os.makedirs(save_folder, exist_ok=True)

# Base URL and compound ID range
base_url = "https://cb.imsc.res.in/imppat/images/3D/SDF/"
start_id = 1
end_id = 17967

# Loop through the compound IDs
for i in range(start_id, end_id + 1):
    # Format the compound ID with leading zeros
    compound_id = f"IMPHY{i:06}"

    # Construct the file URL and save path
    file_url = f"{base_url}{compound_id}_3D.sdf"
    save_path = os.path.join(save_folder, f"{compound_id}_3D.sdf")

    try:
        # Download the file
        response = requests.get(file_url)
        response.raise_for_status()  # Raise an error for HTTP issues

        # Save the file
        with open(save_path, "wb") as file:
            file.write(response.content)

        print(f"Downloaded and saved: {compound_id}")
    except requests.exceptions.RequestException as e:
        print(f"Failed to download {compound_id}: {e}")
