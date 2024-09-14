import pandas as pd
import numpy as np

# Load the CSV files
values_df = pd.read_csv('values.csv')
oil_viscosity_df = pd.read_csv('oil_viscosity.csv')
porosity_df = pd.read_csv('porosity.csv')
permittivity_df = pd.read_csv('permittivity.csv')
tau_porosity_df = pd.read_csv('tau_porosity.csv')
growth_results_df = pd.read_csv('growth_results.csv')  # Load growth_results.csv

# Extract the API value from values.csv
api_value = values_df['API'].iloc[0]  # API is in the first row

# Part 1: API and Viscosity Handling
# Find the nearest API_gravity value to the API in values.csv
oil_viscosity_df['difference'] = np.abs(oil_viscosity_df['API_gravity'] - api_value)
nearest_row = oil_viscosity_df.loc[oil_viscosity_df['difference'].idxmin()]

# Extract the Viscosity_100F value corresponding to the nearest API_gravity
nearest_api_gravity = nearest_row['API_gravity']
viscosity_value = nearest_row['Viscosity_100F']

# Part 2: Porosity Handling
# Check if the Porosity value is N/A
if pd.isna(values_df['Porosity'].iloc[0]):
    # If Porosity is N/A, check the clay content
    clay_content_value = values_df['Clay Content'].iloc[0]
    
    # Find the nearest Polyfit Clay Content (x_smooth) value in Porosity.csv
    porosity_df['difference'] = np.abs(porosity_df['Polyfit Clay Content (x_smooth)'] - clay_content_value)
    nearest_porosity_row = porosity_df.loc[porosity_df['difference'].idxmin()]

    # Extract the Polyfit Porosity (y_smooth) value corresponding to the nearest Polyfit Clay Content (x_smooth)
    nearest_clay_content = nearest_porosity_row['Polyfit Clay Content (x_smooth)']
    porosity_value = nearest_porosity_row['Polyfit Porosity (y_smooth)']
else:
    # If Porosity has a value, use it directly
    porosity_value = values_df['Porosity'].iloc[0]
    nearest_clay_content = 'N/A'  # Not applicable as we don't need to compute

tau_porosity_df['difference'] = np.abs(tau_porosity_df['Porosity'] - porosity_value)
nearest_tau_porosity_row = tau_porosity_df.loc[tau_porosity_df['difference'].idxmin()]
nearest_tau_porosity = nearest_tau_porosity_row['Porosity']
tau_value = nearest_tau_porosity_row['Tortuosity']

# Part 3: Relative Permittivity Handling
# Check if the Relative Permittivity value is N/A
if pd.isna(values_df['Relative Permittivity'].iloc[0]):
    # If Relative Permittivity is N/A, check the Frequency value
    frequency_value = values_df['Frequency'].iloc[0]
    
    # Ensure frequency_value is a float
    try:
        frequency_value = float(frequency_value)
    except ValueError:
        print("Error: Frequency value in values.csv is not a valid float.")
        frequency_value = np.nan

    # Ensure 'Frequency (Hz)' column is numeric
    permittivity_df['Frequency (Hz)'] = pd.to_numeric(permittivity_df['Frequency (Hz)'], errors='coerce')

    # Find the nearest Frequency (Hz) value in Permittivity.csv
    permittivity_df['difference'] = np.abs(permittivity_df['Frequency (Hz)'] - frequency_value)
    nearest_permittivity_row = permittivity_df.loc[permittivity_df['difference'].idxmin()]

    # Extract the Relative Permittivity value corresponding to the nearest Frequency (Hz)
    nearest_frequency = nearest_permittivity_row['Frequency (Hz)']
    relative_permittivity_value = nearest_permittivity_row['Relative Permittivity']
else:
    # If Relative Permittivity has a value, use it directly
    relative_permittivity_value = values_df['Relative Permittivity'].iloc[0]
    nearest_frequency = 'N/A'  # Not applicable as we don't need to compute

# Part 4: Growth Rate (k) Handling
# Extract the Growth Rate (k) value from growth_results.csv
growth_rate_k_value = growth_results_df['Growth Rate (k)'].iloc[0]  # Assuming the Growth Rate (k) is in the first row

# Save the result to data.csv
result_df = pd.DataFrame({
    'API': [api_value], 
    'Nearest_API_gravity': [nearest_api_gravity], 
    'Viscosity_100F': [viscosity_value],
    'Porosity': [porosity_value],
    'Nearest_Clay_Content': [nearest_clay_content],
    'Relative Permittivity': [relative_permittivity_value],
    'Nearest_Frequency': [nearest_frequency],
    'Nearest_Porosity_in_tau': [nearest_tau_porosity],  # Save the nearest porosity value from tau_porosity.csv
    'Tortuosity': [tau_value],  # Save the Tortuosity value
    'Growth Rate (k)': [growth_rate_k_value]  # Save the Growth Rate (k) value from growth_results.csv
})

result_df.to_csv('data.csv', index=False)

print(f"API: {api_value}, Nearest API_gravity: {nearest_api_gravity}, Viscosity_100F: {viscosity_value}, "
      f"Porosity: {porosity_value}, Relative Permittivity: {relative_permittivity_value}, Growth Rate (k): {growth_rate_k_value} saved to data.csv.")
