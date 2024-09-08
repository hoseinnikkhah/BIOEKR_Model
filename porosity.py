import matplotlib.pyplot as plt
import pandas as pd

# Values for the clay content percentage
clay_content_percentage = [
    0, 3.45, 6.9, 10.34, 13.79, 
    17.24, 20.69, 24.14, 27.59, 
    31.03, 34.48, 37.93, 41.38, 
    44.83, 48.28, 51.72, 55.17, 
    58.62, 62.07, 65.52, 68.97, 
    72.41, 75.86, 79.31, 82.76, 
    86.21, 89.66, 93.10, 96.55, 100
]

# Porosity values, initially in percentage
porosity_percentage = [
    25, 24, 23, 22, 21, 
    20, 18, 17, 16, 16, 
    15, 15, 14, 13, 12, 
    11, 17, 20, 24, 25, 
    26, 30, 33, 35, 38, 
    40, 46, 50, 55, 61
]

# Convert porosity percentages to fractions for plotting
porosity_fraction = [value / 100 for value in porosity_percentage]

# Plotting Porosity as a function of Clay Content
plt.figure(figsize=(10, 6))
plt.plot(clay_content_percentage, porosity_fraction, marker='o', color='b', linestyle='-', linewidth=2)
plt.title('Porosity as a Function of Clay Content', fontsize=16)
plt.xlabel('Clay Content (%)', fontsize=14)
plt.ylabel('Porosity (Fraction)', fontsize=14)
plt.grid(True)
plt.savefig('porosity_vs_clay_content.png')
plt.show()

# Save the clay content and porosity data to a CSV file
data_dict = {
    'Clay Content (%)': clay_content_percentage, 
    'Porosity (Fraction)': porosity_fraction
}
data_frame = pd.DataFrame(data_dict)
data_frame.to_csv('porosity_data.csv', index=False)