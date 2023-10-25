# %%
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np

# %%
# Read the text file
with open('L1L2_image.txt', 'r') as file:
    lines = file.readlines()

# %%
# Map characters '.' and 'x' to colors
color_map = {
    '.': [1.0, 1.0, 1.0],  # White color for '.'
    'x': [0.2, 0.2, 0.2]   # Black color for 'x'
}

# Convert characters to colors and create a 2D numpy array
pixels = np.array([[color_map[char] for char in line.strip()] for line in lines])
pixels

# %%
# Display the image using matplotlib
fig, ax = plt.subplots(figsize=(6, 6))
ax.imshow(pixels)
ax.axis('off')  # Turn off axis for better visualization

# Draw a rectangle to indicate the region
ax.axvline(x=24.5, color='y', linewidth=1, linestyle='--')
ax.axhline(y=24.5, color='y', linewidth=1, linestyle='--')
ax.text(12, -3, 'L1', fontsize=16, color='black', ha='center', va='center')
ax.text(37, -3, 'L2', fontsize=16, color='black', ha='center', va='center')
ax.text(-3, 12, 'L1', fontsize=16, color='black', ha='center', va='center', rotation=90)
ax.text(-3, 37, 'L2', fontsize=16, color='black', ha='center', va='center', rotation=90)
# %%
fig.savefig('L1L2_image.pdf', bbox_inches='tight', pad_inches=0.2)

# %%
