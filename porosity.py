import numpy as np
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

# Fit a polynomial of degree 3 (you can adjust the degree)
degree = 3
coefficients = np.polyfit(clay_content_percentage, porosity_fraction, degree)

# Generate a smooth curve by using the polynomial coefficients
polynomial = np.poly1d(coefficients)

# Generate a range of x values for a smooth curve
x_smooth = np.linspace(min(clay_content_percentage), max(clay_content_percentage), 500)
y_smooth = polynomial(x_smooth)

# Plotting Porosity as a function of Clay Content
plt.figure(figsize=(10, 6))
plt.plot(clay_content_percentage, porosity_fraction, 'o', label='Original Data', color='b')
plt.plot(x_smooth, y_smooth, label=f'Polynomial Fit (Degree {degree})', color='r', linestyle='-', linewidth=2)
plt.title('Porosity as a Function of Clay Content (Smoothed)', fontsize=16)
plt.xlabel('Clay Content (%)', fontsize=14)
plt.ylabel('Porosity (Fraction)', fontsize=14)
plt.grid(True)
plt.legend()
plt.savefig('porosity_vs_clay_content_smoothed.png')
plt.show()

max_length = max(len(clay_content_percentage), len(x_smooth))

clay_content_percentage_padded = clay_content_percentage + [np.nan] * (max_length - len(clay_content_percentage))
porosity_fraction_padded = porosity_fraction + [np.nan] * (max_length - len(porosity_fraction))
x_smooth_padded = list(x_smooth) + [np.nan] * (max_length - len(x_smooth))
y_smooth_padded = list(y_smooth) + [np.nan] * (max_length - len(y_smooth))

data_dict = {
    'Clay Content (%)': clay_content_percentage_padded,
    'Porosity (Fraction)': porosity_fraction_padded,
    'Polyfit Clay Content (x_smooth)': x_smooth_padded,
    'Polyfit Porosity (y_smooth)': y_smooth_padded
}

data_frame = pd.DataFrame(data_dict)
data_frame.to_csv('porosity.csv', index=False)
