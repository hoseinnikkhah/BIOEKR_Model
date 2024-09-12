import numpy as np
import matplotlib.pyplot as plt
from numpy.polynomial.polynomial import Polynomial

# Hypothetical data based on the provided plot image
porosity = np.linspace(0.01, 1, 100)  # Porosity from 0 to 1
# Tortuosity decreases as porosity increases (similar to the shape of the curve in the image)
tortuosity = 1 / porosity + 1

# Perform a polynomial fit (degree 2, 3, or higher can be tried based on how well it fits)
degree = 3
coefficients = np.polyfit(porosity, tortuosity, degree)
polynomial = np.poly1d(coefficients)

# Generate fitted values
tortuosity_fit = polynomial(porosity)

# Plotting the original data and the fitted curve
plt.scatter(porosity, tortuosity, label='Original Data')
plt.plot(porosity, tortuosity_fit, color='red', label=f'Poly fit (degree {degree})')
plt.xlabel('Porosity ϕ')
plt.ylabel('Tortuosity τ')
plt.title('Tortuosity vs. Porosity with Polynomial Fit')
plt.legend()
plt.grid(True)
plt.show()
