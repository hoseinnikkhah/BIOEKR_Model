import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

# Generate data
porosity = np.linspace(0.01, 1, 100)
tortuosity = 0.85 / (porosity + 0.1)

# Save the data to a CSV file
data = pd.DataFrame({'Porosity': porosity, 'Tortuosity': tortuosity})
data.to_csv('tau_porosity.csv', index=False)

# Plot the data
plt.plot(porosity, tortuosity, label='Gathered Data', color='black')
plt.xlabel('Porosity ϕ')
plt.ylabel('Tortuosity τ')
plt.title('Tortuosity vs. Porosity')
plt.grid(True)
plt.legend()
plt.show()
