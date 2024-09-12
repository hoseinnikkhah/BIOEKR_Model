import pandas as pd
import numpy as np

# Load the CSV files
values_df = pd.read_csv('values.csv')
oil_viscosity_df = pd.read_csv('oil_viscosity.csv')

# Extract the API value from values.csv
api_value = values_df['API'].iloc[0]  # Assuming API is in the first row

# Find the nearest API_gravity value to the API in values.csv
oil_viscosity_df['difference'] = np.abs(oil_viscosity_df['API_gravity'] - api_value)
nearest_row = oil_viscosity_df.loc[oil_viscosity_df['difference'].idxmin()]

# Extract the Viscosity_100F value corresponding to the nearest API_gravity
nearest_api_gravity = nearest_row['API_gravity']
viscosity_value = nearest_row['Viscosity_100F']

# Save the result to data.csv
result_df = pd.DataFrame({'API': [api_value], 'Nearest_API_gravity': [nearest_api_gravity], 'Viscosity_100F': [viscosity_value]})
result_df.to_csv('data.csv', index=False)

print(f"API: {api_value}, Nearest API_gravity: {nearest_api_gravity}, Viscosity_100F: {viscosity_value} saved to data.csv.")
