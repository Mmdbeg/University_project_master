import georinex as gr
import pandas as pd
import numpy as np 
import re
import json  
import os 
from datetime import datetime
from scipy.interpolate import lagrange


# 1. Define the calculation function (using georinex column names)
def calculate_ecef(row, mu, omega_e):
    try:
        # EXACT mapping from your column list
        sqA = row['sqrtA']
        e   = row['Eccentricity']  # Fixed
        M0  = row['M0']
        dn  = row['DeltaN']        # Fixed
        toe = row['Toe']
        omg = row['omega']
        Om0 = row['Omega0']
        i0  = row['Io']            # Fixed
        idot = row['IDOT']
        OmDot = row['OmegaDot']
        
        # Corrections
        cus, cuc = row['Cus'], row['Cuc']
        cis, cic = row['Cis'], row['Cic']
        crs, crc = row['Crc'], row['Crs']

        # 1. Time calculation
        t = row['time'].timestamp() % 604800 
        tk = t - toe
        
        if tk > 302400: tk -= 604800
        if tk < -302400: tk += 604800

        # 2. Orbit math
        A = sqA**2
        n = np.sqrt(mu/(A**3)) + dn
        Mk = M0 + n * tk
        
        # Kepler's Equation
        Ek = Mk
        for _ in range(10):
            Ek = Ek - (Ek - e * np.sin(Ek) - Mk) / (1 - e * np.cos(Ek))
            
        nu_k = np.arctan2(np.sqrt(1 - e**2) * np.sin(Ek), np.cos(Ek) - e)
        phi_k = nu_k + omg
        
        # 3. Perturbations
        uk = phi_k + cus*np.sin(2*phi_k) + cuc*np.cos(2*phi_k)
        rk = A*(1 - e*np.cos(Ek)) + crc*np.cos(2*phi_k) + crs*np.sin(2*phi_k)
        ik = i0 + idot*tk + cic*np.cos(2*phi_k) + cis*np.sin(2*phi_k)
        
        # 4. Orbital plane to ECEF
        x_p = rk * np.cos(uk)
        y_p = rk * np.sin(uk)
        
        Omk = Om0 + (OmDot - omega_e)*tk - omega_e*toe
        
        X = x_p * np.cos(Omk) - y_p * np.cos(ik) * np.sin(Omk)
        Y = x_p * np.sin(Omk) + y_p * np.cos(ik) * np.cos(Omk)
        Z = y_p * np.sin(ik)
        
        return pd.Series([X, Y, Z], index=['X', 'Y', 'Z'])

    except Exception as err:
        return pd.Series([np.nan, np.nan, np.nan], index=['X', 'Y', 'Z'])




# --- Constants & Configuration ---

PARAM_KEYS = [
    ["af0", "af1", "af2"],
    ["IODE_IODnav", "Crs", "Delta_n", "M0"],
    ["Cuc", "e", "Cus", "sqrtA"],
    ["Toe", "Cic", "Omega0", "Cis"],
    ["i0", "Crc", "omega", "Omega_dot"],
    ["IDOT", "Week", "SISA_Acc", "Health"],
    ["TGD_BGDe5a", "IODC_BGDe5b", "Toc", "Fit"],
    ["Spare1", "Spare2"]
]

# --- Helper Functions ---

def extract_all_numbers_from_line(line, is_line1=False):
    pattern = r'[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?'
    numeric_strings = re.findall(pattern, line.replace('D', 'E')) # Handle D notation
    numbers = [float(p) for p in numeric_strings]
    
    if is_line1:
        return numbers[-3:] if len(numbers) >= 3 else []
    return numbers

def parse_nav_block(data_block):
    lines = data_block.split('\n')
    if len(lines) < 8: return None
    
    # --- 1. Extract SV_ID and Time (Line 1) ---
    header = lines[0]
    sv_id = header[0:3].strip() # Extracts 'E02'
    
    try:
        # Based on your structure: E02 YYYY MM DD HH mm SS
        # We pick specific indices to avoid the '02' in 'E02'
        year   = int(header[4:8])   # 2025
        month  = int(header[9:11])  # 04
        day    = int(header[12:14]) # 10
        hour   = int(header[15:17]) # 00
        minute = int(header[18:20]) # 00
        second = int(header[21:23]) # 00
        
        dt = datetime(year, month, day, hour, minute, second)
    except (ValueError, IndexError):
        # Fallback: find all numbers and skip the first one (which is the '02' from E02)
        nums = re.findall(r'\d+', header)
        if len(nums) >= 7:
            dt = datetime(int(nums[1]), int(nums[2]), int(nums[3]), 
                          int(nums[4]), int(nums[5]), int(nums[6]))
        else:
            dt = None

    # --- 2. Extract Numerical Parameters (Scientific Notation) ---
    # Galileo often has no spaces between negative numbers: e.g. 0.000e+00-6.519e-08
    # This regex splits them correctly
    def get_numbers(line):
        return [float(x) for x in re.findall(r'[+-]?\d+(?:\.\d+)?(?:[eEdD][+-]?\d+)?', line.replace('D', 'E'))]

    # Line 1 Clock Params (the 3 numbers after the date)
    clock_params = get_numbers(header[23:])
    
    # Lines 2-8 Orbital Params
    body_params = []
    for line in lines[1:8]:
        body_params.extend(get_numbers(line))
        
    # --- 3. Build the Record ---
    record = {
        'SV_ID': sv_id, 
        'System': 'Galileo',
        'time': dt 
    }
    
    # Map Clock Params (af0, af1, af2)
    for i, val in enumerate(clock_params):
        if i < 3: record[PARAM_KEYS[0][i]] = val
        
    # Map Orbital Params
    for i, val in enumerate(body_params):
        if i < 28:
            row_idx = (i // 4) + 1
            col_idx = i % 4
            record[PARAM_KEYS[row_idx][col_idx]] = val
            
    return record

def process_rinex_to_df(rinex_content):
    header_end = rinex_content.find("END OF HEADER")
    if header_end == -1: return pd.DataFrame()

    data_section = rinex_content[header_end + len("END OF HEADER"):]
    # Split by satellite identifier (Gxx, Exx, etc.) at the start of a line
    block_pattern = r"^\s*([GE]\d{2}\s\d{4}\s)"
    parts = re.split(block_pattern, data_section, flags=re.MULTILINE)
    
    records = []
    for i in range(1, len(parts), 2):
        if i + 1 < len(parts):
            full_block = parts[i] + parts[i+1]
            data = parse_nav_block(full_block)
            if data: records.append(data)
            
    df = pd.DataFrame(records)
    # Add calculated Semi-major axis
    if 'sqrtA' in df.columns:
        df['A'] = df['sqrtA']**2
    return df



# --- Kepler Solver (Your Logic) ---

def solve_kepler(row):
    M0 = row['M0']
    e = row['e']
    # Start E with M0
    Ek = M0
    for _ in range(10):
        Ek = Ek - (Ek - e * np.sin(Ek) - M0) / (1 - e * np.cos(Ek))
    return Ek

def calculate_galileo_ecef(row):
    try:
        # --- Constants (Galileo GTRF) ---
        mu = 3.986004418e14
        omega_e = 7.2921151467e-5
        
        # --- Extract values with fallbacks ---
        # Using .get() prevents the function from crashing if a column is missing
        toe = row.get('Toe')
        A = row.get('A')
        e = row.get('e')
        M0 = row.get('M0')
        dn = row.get('Delta_n')
        omg = row.get('omega')
        Om0 = row.get('Omega0')
        OmDot = row.get('Omega_dot')
        i0 = row.get('i0')
        idot = row.get('IDOT')

        # --- 1. Time Calculation ---
        t = row['time'].timestamp() % 604800 
        tk = t - toe
        if tk > 302400: tk -= 604800
        if tk < -302400: tk += 604800

        # --- 2. Mean Motion & Mean Anomaly ---
        n = np.sqrt(mu / (A**3)) + dn
        Mk = M0 + n * tk

        # --- 3. Solve Kepler ---
        Ek = Mk
        for _ in range(10):
            Ek = Ek - (Ek - e * np.sin(Ek) - Mk) / (1 - e * np.cos(Ek))

        # --- 4. Argument of Latitude ---
        nu = np.arctan2(np.sqrt(1 - e**2) * np.sin(Ek), np.cos(Ek) - e)
        phi = nu + omg

        # --- 5. Perturbations ---
        # Ensure these match your PARAM_KEYS names exactly!
        u = phi + row['Cus']*np.sin(2*phi) + row['Cuc']*np.cos(2*phi)
        r = A*(1 - e*np.cos(Ek)) + row['Crs']*np.sin(2*phi) + row['Crc']*np.cos(2*phi)
        i = i0 + idot*tk + row['Cis']*np.sin(2*phi) + row['Cic']*np.cos(2*phi)

        # --- 6. Orbital Plane to ECEF ---
        x_p = r * np.cos(u)
        y_p = r * np.sin(u)
        Omk = Om0 + (OmDot - omega_e)*tk - (omega_e * toe)

        X = x_p * np.cos(Omk) - y_p * np.cos(i) * np.sin(Omk)
        Y = x_p * np.sin(Omk) + y_p * np.cos(i) * np.cos(Omk)
        Z = y_p * np.sin(i)

        return pd.Series({'X': X, 'Y': Y, 'Z': Z})

    except Exception as err:
        # This will tell you EXACTLY what is missing
        print(f"Error in satellite {row.get('SV_ID')}: {err}")
        return pd.Series({'X': np.nan, 'Y': np.nan, 'Z': np.nan})