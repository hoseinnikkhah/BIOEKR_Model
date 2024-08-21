import numpy as np
import matplotlib.pyplot as plt

# Time data ranging from 1 to 35 days
time_data = np.arange(1, 36)

# Concentration values as an array
concentration_data = np.array([
    1.64872127070013, 1.68301736879343, 1.71802688180127, 1.75376464992030, 1.79024582204753,
    1.82748586220180, 1.86550055607883, 1.90430601774260, 1.94391869645598, 1.98435538365342,
    2.02563322005864, 2.06776970295043, 2.11078269357961, 2.15469042474020, 2.19951150849817,
    2.24526494408088, 2.29197012593066, 2.33964685192599, 2.38831533177357, 2.43799619557505,
    2.48871050257191, 2.54047975007231, 2.59332588256355, 2.64727130101414, 2.70233887236937,
    2.75855193924439, 2.81593432981889, 2.87451036793768, 2.93430488342127, 2.99534322259109,
    3.05765125901345, 3.12125540446716, 3.18618262013922, 3.25246042805340, 3.32011692273655
])

# Compute the natural log of the concentration data
ln_concentration = np.log(concentration_data)

# Create a plot of the logarithmic concentration over time
plt.figure(figsize=(8, 5))
plt.plot(time_data, ln_concentration, color='purple', marker='o', linestyle='--', label='Log(Concentration)')
plt.xlabel('Time (days)')
plt.ylabel('ln(Microbial Concentration) (mg/L)')
plt.title('Logarithmic Microbial Growth Over Time')
plt.grid(True, which='both', linestyle=':')
plt.legend()
plt.tight_layout()
plt.show()

# Perform linear regression and obtain slope and intercept
linear_coefficients = np.polyfit(time_data, ln_concentration, 1)
growth_rate = linear_coefficients[0]
y_intercept = linear_coefficients[1]

# Display the results
print(f"Growth rate (k): {growth_rate:.4f}")
print(f"Y-Intercept: {y_intercept:.4f}")
