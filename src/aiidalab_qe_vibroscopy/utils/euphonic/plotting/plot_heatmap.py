# to be able to plot:
# pip install pandas json matplotlib


import pandas as pd
import json
import matplotlib.pyplot as plt

# Load the heatmap data from a CSV file
csv_file = "INS_structure_factor_22.csv"  # Change this to your CSV file path
df = pd.read_csv(csv_file, index_col=0)

json_file = "INS_metadata_22.json"  # Change this to your JSON file path
with open(json_file, "r") as f:
    metadata = json.load(f)

# Create the heatmap
plt.figure(figsize=(10, 8))
plt.imshow(df.values, cmap="cividis", aspect="auto")

plt.ylabel = metadata["ylabel"]

# Add color bar
plt.colorbar()

# Add x and y axis labels
if "ticks_labels" in metadata:
    plt.xticks(metadata["ticks_positions"], metadata["ticks_labels"])

if "spectrum_type" in metadata:
    plt.title(f"Inelastic neutron scattering data - {metadata['spectrum_type']}")
else:
    plt.title("Inelastic neutron scattering data")

plt.gca().invert_yaxis()
# Show the plot
plt.tight_layout()
plt.show()
