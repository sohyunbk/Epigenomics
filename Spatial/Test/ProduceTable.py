import os

# Debugging statements
print("Current working directory:", os.getcwd())
output_path = "/scratch/sb14489/Test.txt"
outfile = open(output_path, "w")
outfile.write("Test")
outfile.close()

# Check if the file was created successfully
if os.path.exists(output_path):
    print("File created successfully:", output_path)
else:
    print("Failed to create the file:", output_path)
