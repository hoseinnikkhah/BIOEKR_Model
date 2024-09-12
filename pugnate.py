import subprocess

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
    tortuosity = None
    has_porosity = input("Do we have a porosity value? (yes/no): ").strip().lower()

    if has_porosity == 'yes':
        porosity = input("What is the porosity value? (Enter a number): ")
    else:
        tortuosity = input("What is the tortuosity value? (Enter a number): ")

    # 3. Ask if the crude oil relative permittivity is determined
    relative_permittivity = None
    frequency = None
    has_relative_permittivity = input("Is crude oil relative permittivity determined? (yes/no): ").strip().lower()

    if has_relative_permittivity == 'yes':
        relative_permittivity = input("What is the relative permittivity value? (Enter a number): ")
    else:
        frequency = input("What is the frequency? (Enter a value between 10^2 and 10^10 Hz): ")

    # Return all gathered inputs
    return API, porosity, tortuosity, relative_permittivity, frequency

# Main function
def main():
    # Run the Python scripts
    run_scripts()

    # Gather inputs from the user
    API, porosity, tortuosity, relative_permittivity, frequency = gather_inputs()

    # Output the gathered data for further use
    print(f"\nCollected Data:")
    print(f"API: {API}")
    print(f"Porosity: {porosity if porosity else 'Not provided'}")
    print(f"Tortuosity: {tortuosity if tortuosity else 'Not provided'}")
    print(f"Relative Permittivity: {relative_permittivity if relative_permittivity else 'Not provided'}")
    print(f"Frequency: {frequency if frequency else 'Not provided'}")

if __name__ == "__main__":
    main()
