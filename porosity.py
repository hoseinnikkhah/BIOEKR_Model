import matplotlib.pyplot as plt
import pandas as pd

# Clay content values
clay_content = [0, 3.44827586206897, 6.89655172413793, 10.3448275862069, 13.7931034482759, 
                17.2413793103448, 20.6896551724138, 24.1379310344828, 27.5862068965517, 
                31.0344827586207, 34.4827586206897, 37.9310344827586, 41.3793103448276, 
                44.8275862068966, 48.2758620689655, 51.7241379310345, 55.1724137931034, 
                58.6206896551724, 62.0689655172414, 65.5172413793104, 68.9655172413793, 
                72.4137931034483, 75.8620689655172, 79.3103448275862, 82.7586206896552, 
                86.2068965517241, 89.6551724137931, 93.1034482758621, 96.5517241379310, 100]

# Porosity values
porosity = [25, 24, 23, 22, 21, 20, 18, 17, 16, 16, 15, 15, 14, 13, 12, 11, 17, 20, 
            24, 25, 26, 30, 33, 35, 38, 40, 46, 50, 55, 61]

# Plotting Porosity vs. Clay Content
plt.figure(figsize=(10, 6))
plt.plot(clay_content, porosity, marker='o')
plt.title('Porosity vs. Clay Content')
plt.xlabel('Clay Content (%)')
plt.ylabel('Porosity')
plt.grid(True)
plt.savefig('porosity_vs_clay_content.png')
plt.show()

# Saving the data to CSV
data = {'Clay Content (%)': clay_content, 'n (Porosity)': porosity}
df = pd.DataFrame(data)
df.to_csv('porosity_set.csv', index=False)
