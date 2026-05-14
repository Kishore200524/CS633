import matplotlib.pyplot as plt
import numpy as np
import csv

data_small = []
data_large = []

filename = 'timing_data.txt'

try:
    with open(filename, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            # Format: [P, Size, T1, T2, T3, T4, T5]
            p = int(row[0].strip())
            size = row[1].strip().lower()
            times = [float(x) for x in row[2:]]
            
            if size == 'small':
                data_small.append(times)
            elif size == 'large':
                data_large.append(times)
                
    print(f"Successfully read data for {len(data_small)} process configurations.")

except FileNotFoundError:
    print(f"Error: Could not find '{filename}'. Make sure it is in the same folder.")
    exit()


fig, ax = plt.subplots(figsize=(10, 6))

positions = np.array([1, 2, 3])
width = 0.3

# Plot Small M 
bp1 = ax.boxplot(data_small, positions=positions - 0.2, widths=width, patch_artist=True,
                 boxprops=dict(facecolor="skyblue", color="blue"),
                 medianprops=dict(color="black"))

# Plot Large M 
bp2 = ax.boxplot(data_large, positions=positions + 0.2, widths=width, patch_artist=True,
                 boxprops=dict(facecolor="orange", color="darkred"),
                 medianprops=dict(color="black"))

# Format
ax.set_title('MPI Performance: Execution Time vs. Processes', fontsize=14)
ax.set_xlabel('Number of Processes (P)', fontsize=12)
ax.set_ylabel('Execution Time (seconds)', fontsize=12)
ax.set_xticks(positions)
ax.set_xticklabels(['8', '16', '32'])
ax.grid(True, linestyle='--', alpha=0.5)

# Legend
ax.legend([bp1["boxes"][0], bp2["boxes"][0]], 
          ['M = 262,144', 'M = 1,048,576'], 
          loc='upper left')

# Save
plt.tight_layout()
output_file = 'performance_boxplot.png'
plt.savefig(output_file, dpi=300)
print(f"Graph saved as '{output_file}'")
