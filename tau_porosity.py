import numpy as np
import matplotlib.pyplot as plt

# Hypothetical data resembling your plot
porosity = np.linspace(0.01, 1, 100)  # Porosity from 0 to 1
# Tortuosity decreases with increasing porosity (with a range up to 8, similar to your plot)
tortuosity = 0.85 / (porosity + 0.2)  # A decreasing curve from 8 to near 1, adjusting for visual similarity

# Plotting the original data as a line
plt.plot(porosity, tortuosity, label='Original Data', color='black')
plt.xlabel('Porosity ϕ')
plt.ylabel('Tortuosity τ')
plt.title('Tortuosity vs. Porosity')
plt.grid(True)
plt.legend()
plt.show()
