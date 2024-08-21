import numpy as np
import matplotlib.pyplot as plt
import csv
import os

# Simulated time data spanning 35 days
days = np.arange(1, 36)

# Given microbial concentration data (mg/L) over the time period
microbial_concentration = np.array([
    1.64872127070013, 1.68301736879343, 1.71802688180127, 1.75376464992030, 1.79024582204753,
    1.82748586220180, 1.86550055607883, 1.90430601774260, 1.94391869645598, 1.98435538365342,
    2.02563322005864, 2.06776970295043, 2.11078269357961, 2.15469042474020, 2.19951150849817,
    2.24526494408088, 2.29197012593066, 2.33964685192599, 2.38831533177357, 2.43799619557505,
    2.48871050257191, 2.54047975007231, 2.59332588256355, 2.64727130101414, 2.70233887236937,
    2.75855193924439, 2.81593432981889, 2.87451036793768, 2.93430488342127, 2.99534322259109,
    3.05765125901345, 3.12125540446716, 3.18618262013922, 3.25246042805340, 3.32011692273655
])

# Calculate the natural logarithm of the concentration data
log_concentration = np.log(microbial_concentration)

# Plotting the log-transformed concentration data over time
plt.figure(figsize=(8, 5))
plt.plot(days, log_concentration, color='purple', marker='o', linestyle='--', label='ln(Concentration)')
plt.xlabel('Time (days)')
plt.ylabel('ln(Microbial Concentration) (mg/L)')
plt.title('Logarithmic Growth of Microbial Concentration')
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()
plt.show()

# Perform linear regression using numpy's polyfit function
linear_fit = np.polyfit(days, log_concentration, 1)
growth_rate_k = linear_fit[0]
y_intercept = linear_fit[1]

# Validate if the file already exists
output_file = 'growth_results.csv'
if os.path.exists(output_file):
    print(f"Warning: '{output_file}' already exists and will be overwritten.")

# Save the growth rate and intercept to a CSV file
with open(output_file, mode='w', newline='') as file:
    csv_writer = csv.writer(file)
    csv_writer.writerow(['Growth Rate (k)', 'Y-Intercept'])
    csv_writer.writerow([growth_rate_k, y_intercept])

# Display calculated growth rate and intercept
print(f"Calculated Growth Rate (k): {growth_rate_k:.4f}")
print(f"Calculated Y-Intercept: {y_intercept:.4f}")
