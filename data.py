import pandas as pd

# Load the CSV files
values_df = pd.read_csv('values.csv')
oil_viscosity_df = pd.read_csv('oil_viscosity.csv')

# Extract the API value from values.csv
api_value = values_df['API'].iloc[0]  # Assuming API is in the first row

# Find the matching Viscosity_100F for the corresponding API_gravity
matching_row = oil_viscosity_df[oil_viscosity_df['API_gravity'] == api_value]

# Check if a match was found
if not matching_row.empty:
    viscosity_value = matching_row['Viscosity_100F'].values[0]
    
    # Save the result to data.csv
    result_df = pd.DataFrame({'API': [api_value], 'Viscosity_100F': [viscosity_value]})
    result_df.to_csv('data.csv', index=False)
    
    print(f"API: {api_value}, Viscosity_100F: {viscosity_value} saved to data.csv.")
else:
    print(f"No matching API_gravity found for API: {api_value}.")
