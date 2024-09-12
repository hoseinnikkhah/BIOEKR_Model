import pandas as pd
import numpy as np

# Load the CSV files
values_df = pd.read_csv('values.csv')
oil_viscosity_df = pd.read_csv('oil_viscosity.csv')
porosity_df = pd.read_csv('porosity.csv')

# Extract the API value from values.csv
api_value = values_df['API'].iloc[0]  # Assuming API is in the first row

# Part 1: API and Viscosity Handling
# Find the nearest API_gravity value to the API in values.csv
oil_viscosity_df['difference'] = np.abs(oil_viscosity_df['API_gravity'] - api_value)
nearest_row = oil_viscosity_df.loc[oil_viscosity_df['difference'].idxmin()]

# Extract the Viscosity_100F value corresponding to the nearest API_gravity
nearest_api_gravity = nearest_row['API_gravity']
viscosity_value = nearest_row['Viscosity_100F']

# Part 2: Porosity Handling
# Check if the porosity value is N/A
if pd.isna(values_df['porosity'].iloc[0]):
    # If porosity is N/A, check the clay content
    clay_content_value = values_df['clay content'].iloc[0]
    
    # Find the nearest Clay Content (%) value in porosity.csv
    porosity_df['difference'] = np.abs(porosity_df['Clay Content (%)'] - clay_content_value)
    nearest_porosity_row = porosity_df.loc[porosity_df['difference'].idxmin()]

    # Extract the Porosity (Fraction) value corresponding to the nearest Clay Content (%)
    nearest_clay_content = nearest_porosity_row['Clay Content (%)']
    porosity_value = nearest_porosity_row['Porosity (Fraction)']
else:
    # If porosity has a value, use it directly
    porosity_value = values_df['porosity'].iloc[0]
    nearest_clay_content = 'N/A'  # Not applicable as we don't need to compute

# Save the result to data.csv
result_df = pd.DataFrame({
    'API': [api_value], 
    'Nearest_API_gravity': [nearest_api_gravity], 
    'Viscosity_100F': [viscosity_value],
    'Porosity': [porosity_value],
    'Nearest_Clay_Content': [nearest_clay_content]
})
result_df.to_csv('data.csv', index=False)

print(f"API: {api_value}, Nearest API_gravity: {nearest_api_gravity}, Viscosity_100F: {viscosity_value}, Porosity: {porosity_value} saved to data.csv.")
