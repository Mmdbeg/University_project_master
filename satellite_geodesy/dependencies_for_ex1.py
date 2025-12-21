import pandas as pd
import numpy as np 
import re
import json  
import os # Added for file path checks

  

# --- Data Structure Mapping and Keys (Fixed Definitions) ---

nav_data_structure = {
    "SV_ID": None,
    "Line1_Clock": {"af0": None, "af1": None, "af2": None},
    "Line2": {"IODE": None, "Crs": None, "Delta_n": None, "M0": None},
    "Line3": {"Cuc": None, "e": None, "Cus": None, "sqrtA": None},
    "Line4": {"Toe": None, "Cic": None, "Omega0": None, "Cis": None},
    "Line5": {"i0": None, "Crc": None, "omega": None, "Omega_dot": None},
    "Line6": {"IDOT": None, "GPS_Week": None, "Accuracy": None, "Health": None},
    "Line7": {"TGD": None, "IODC": None, "Toc": None, "Fit_Interval": None}
}

param_keys = [
    ["af0", "af1", "af2"],
    ["IODE", "Crs", "Delta_n", "M0"],
    ["Cuc", "e", "Cus", "sqrtA"],
    ["Toe", "Cic", "Omega0", "Cis"],
    ["i0", "Crc", "omega", "Omega_dot"],
    ["IDOT", "GPS_Week", "Accuracy", "Health"],
    ["TGD", "IODC", "Toc", "Fit_Interval"],
    ["Reserved1", "Reserved2"]
]

# --- Parsing Utility Functions ---

def extract_all_numbers_from_line(line, is_line1=False):
    """
    Extracts all numbers (including scientific notation and handling no-space 
    negative sign issues) using Regular Expressions.
    """
    pattern = r'[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?'
    numeric_strings = re.findall(pattern, line)
    numbers = [float(p) for p in numeric_strings]
    
    if is_line1:
        if len(numbers) >= 3:
            return numbers[-3:]
        return []
    else:
        return numbers

def extract_keplerian_parameters_regex(data_block):
    """
    Parses a single 8-line navigation data block into a structured dictionary.
    """
    result = nav_data_structure.copy()
    lines = data_block.strip().split('\n')
    all_numeric_fields = []
    
    if not lines or len(lines) < 8:  # this line make out function specific for rinex nav files 
        return None
        
    # Line 1: SV_ID and Clock Parameters
    line1 = lines[0]
    result["SV_ID"] = line1[0:3].strip() # SV_ID 
    
    clock_params = extract_all_numbers_from_line(line1, is_line1=True) # clock
    if len(clock_params) == 3:
        result["Line1_Clock"]["af0"] = clock_params[0]
        result["Line1_Clock"]["af1"] = clock_params[1]
        result["Line1_Clock"]["af2"] = clock_params[2]
    
    # Lines 2 through 8: Orbital Parameters
    for line in lines[1:]:
        all_numeric_fields.extend(extract_all_numbers_from_line(line, is_line1=False))
                
    # Populate the Dictionaries Sequentially

    field_index = 0  # An external counter 

    for i in range(1, len(param_keys)): 
        line_key = f"Line{i+1}"
        current_param_names = param_keys[i]
        
        current_dict = result.get(line_key)
        
        if current_dict is None and i == 7:
            current_dict = result[f"Line{i+1}_Reserved"] = {}
        elif current_dict is None:
            continue

        for param_name in current_param_names:
            if field_index < len(all_numeric_fields):
                current_dict[param_name] = all_numeric_fields[field_index]
                field_index += 1
            else:
                break
    
    return result

def create_keplerian_record(extracted_data):
    """
    Extracts only the core Keplerian elements and rates into a flat dictionary
    suitable for a DataFrame row.
    """
    if not extracted_data:
        return None

    line2 = extracted_data.get('Line2', {})
    line3 = extracted_data.get('Line3', {})
    line4 = extracted_data.get('Line4', {})
    line5 = extracted_data.get('Line5', {})
    line6 = extracted_data.get('Line6', {})
    
    sqrtA = line3.get('sqrtA') # extract sqrtA indivigually 
    
    # Create the single row record
    record = {
        'SV_ID': extracted_data.get('SV_ID', 'N/A'),
        'Epoch_Time': f"{line4.get('Toe')} (sec)", 
        
        # --- The Six Keplerian Elements ---
        'M0 (rad)': line2.get('M0'),
        'e': line3.get('e'),
        'A (m)': sqrtA**2 if sqrtA is not None else None,
        'i0 (rad)': line5.get('i0'),
        'Omega0 (rad)': line4.get('Omega0'),
        'omega (rad)': line5.get('omega'),

        # --- Related Rates and Corrections ---
        'Delta_n (rad/s)': line2.get('Delta_n'),
        'Omega_dot (rad/s)': line5.get('Omega_dot'),
        'IDOT (rad/s)': line6.get('IDOT'),
    }
    return record

# --- Main File Processing Function ---

def process_rinex_file_to_dataframe(rinex_content):
    """
    Processes the entire RINEX navigation file content string, extracts all 
    satellite blocks iteratively, and compiles a single DataFrame.
    """
    
    header_end = rinex_content.find("END OF HEADER")
    if header_end == -1:  # this condition make sure that has "END OF HEADER " 
        print(" can't find 'END OF HEADER' ")
        return pd.DataFrame()

    data_section = rinex_content[header_end + len("END OF HEADER"):]
    
    # Pattern to find the start of a new satellite block (e.g., 'G01 2025')
    block_pattern = r"^\s*([GRECJI]\d{2}\s\d{4}\s\d{2}\s\d{2}\s\d{2}\s\d{2}\s\d{2}\s*)"
    
    blocks = re.split(block_pattern, data_section, flags=re.MULTILINE)
    
    satellite_blocks = []
    
    if len(blocks) > 2:
        for i in range(1, len(blocks), 2):
            if i + 1 < len(blocks):
                full_block = blocks[i] + blocks[i+1]
                satellite_blocks.append(full_block.strip())
            
    all_records = []
    
    print(f"{len(satellite_blocks)} satellite found ")

# ------------------------------------------------------------------

    for block in satellite_blocks:
        parsed_data = extract_keplerian_parameters_regex(block)
        
        if parsed_data:
            keplerian_record = create_keplerian_record(parsed_data)
            
            if keplerian_record:
                all_records.append(keplerian_record)
                
    # Final DataFrame Creation
    df_kaplerian = pd.DataFrame(all_records)
    if not df_kaplerian.empty:
        df_kaplerian = df_kaplerian.set_index('SV_ID')

    return df_kaplerian

def newton_raphson(M0, e, E0, tolerance=1e-12):
    En = E0
    for i in range(100):
        # f(E) = E - e*sin(E) - M
        # f'(E) = 1 - e*cos(E)
        # Newton-Raphson: En_new = En - f(En)/f'(En)
        En_new = En - ((En - e * np.sin(En) - M0) / (1 - e * np.cos(En)))

        # Check the difference between THIS step and the PREVIOUS step
        if abs(En_new - En) <= tolerance:
            return En_new
            
        En = En_new  # Update En for the next iteration

    return En  # Return the best result found if it doesn't converge in 100 steps
   
def series_eccentric_anomaly(M, e):
    """
    Calculates Eccentric Anomaly (E) using the series expansion formula

    """
    # Pre-calculating powers of e for efficiency
    e2, e3, e4, e5, e6, e7 = e**2, e**3, e**4, e**5, e**6, e**7
    
    # Calculating each trigonometric term based on the image
    term1 = (e - (1/8)*e3 + (1/192)*e5 - (1/9216)*e7) * np.sin(M)
    term2 = (0.5*e2 - (1/6)*e4 + (1/48)*e6) * np.sin(2*M)
    term3 = ((3/8)*e3 - (27/128)*e5 + (243/5120)*e7) * np.sin(3*M)
    term4 = ((1/3)*e4 - (4/15)*e6) * np.sin(4*M)
    term5 = ((125/384)*e5 - (3125/9216)*e7) * np.sin(5*M)
    term6 = (27/80)*e6 * np.sin(6*M)
    term7 = -(16807/46080)*e7 * np.sin(7*M)
    
    # Summing it all up: E = M + terms
    E = M + term1 + term2 + term3 + term4 + term5 + term6 + term7
    return E



