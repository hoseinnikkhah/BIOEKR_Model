import subprocess
import csv

# Run the specified Python files
def run_scripts():
    scripts = ['oil_viscosity.py', 'permittivity.py', 'porosity.py', 'tau_porosity.py', 'microbial_growth.py']
    for script in scripts:
        try:
            subprocess.run(['python3', script], check=True)
            print(f"Successfully ran {script}")
        except subprocess.CalledProcessError as e:
            print(f"Error running {script}: {e}")

# Gather input from the user
def gather_inputs():
    # 1. Ask for the API value of the crude oil
    API = input("What API is the crude oil? (Enter a number): ")

    # 2. Ask if we have a porosity value
    porosity = None
    clay_content = None
    has_porosity = input("Do we have a porosity value? (yes/no): ").strip().lower()

    if has_porosity == 'yes':
        porosity = input("What is the porosity value? (Enter a number): ")
    else:
        clay = input("What is the clay_content value? (Enter a number between 0 and 1): ")

    # 3. Ask if the crude oil relative permittivity is determined
    relative_permittivity = None
    frequency = None
    has_relative_permittivity = input("Is crude oil relative permittivity determined? (yes/no): ").strip().lower()

    if has_relative_permittivity == 'yes':
        relative_permittivity = input("What is the relative permittivity value? (Enter a number): ")
    else:
        frequency = input("What is the frequency? (Enter a value between 10^2 and 10^10 Hz): ")

    # Return all gathered inputs
    return API, porosity, clay_content, relative_permittivity, frequency

# Save data to a CSV file
def save_to_csv(data):
    filename = 'values.csv'
    fieldnames = ['API', 'Porosity', 'Tortuosity', 'Relative Permittivity', 'Frequency']

    # Open the CSV file and write the data
    with open(filename, mode='w', newline='') as file:
        writer = csv.DictWriter(file, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerow({
            'API': data['API'],
            'Porosity': data['Porosity'] if data['Porosity'] else 'N/A',
            'Tortuosity': data['Tortuosity'] if data['Tortuosity'] else 'N/A',
            'Relative Permittivity': data['Relative Permittivity'] if data['Relative Permittivity'] else 'N/A',
            'Frequency': data['Frequency'] if data['Frequency'] else 'N/A'
        })
    print(f"Data saved to {filename}")

# Main function
def main():
    # Run the Python scripts
    run_scripts()

    # Gather inputs from the user
    API, porosity, clay_content, relative_permittivity, frequency = gather_inputs()

    # Store collected data in a dictionary
    data = {
        'API': API,
        'Porosity': porosity,
        'Tortuosity': clay_content,
        'Relative Permittivity': relative_permittivity,
        'Frequency': frequency
    }

    # Output the gathered data
    print(f"\nCollected Data:")
    for key, value in data.items():
        print(f"{key}: {value if value else 'Not provided'}")

    # Save data to CSV
    save_to_csv(data)

if __name__ == "__main__":
    main()
