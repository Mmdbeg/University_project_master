from dependencies_for_ex3 import *


# LOADING FILES PATH 
file_path_GPS = "GPS/GPS_only_nav.rnx"
file_path_galileo = "galileo/GAL_only_nav.rnx"
file_path_beidu = "beidu/BDS_only_nav.rnx"


# ... (your loading code GPS)
Gps_nav = gr.load(file_path_GPS)
df_gps_only = Gps_nav.to_dataframe().reset_index()




# #    ... (your loading code Beidu)
# Beidu_nav = gr.load(file_path_beidu)
# df_bds_only = Beidu_nav.to_dataframe().reset_index()
                                                   

#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
  #                                                            wgs-84 *** GPS ellipsoid
ecef_results_gps = df_gps_only.apply(calculate_ecef, axis=1, args=(3.986005e14, 7.2921151467e-5))
df_gps_only = pd.concat([df_gps_only, ecef_results_gps], axis=1)

#Remove all rows where X, Y, or Z could not be calculated
df_clean = df_gps_only.dropna(subset=['X', 'Y', 'Z'])
print(f"Original rows: {len(df_gps_only)}")
print(f"Cleaned rows: {len(df_clean)}")
print(df_clean[['sv', 'time', 'X', 'Y', 'Z']])



#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
# #                                                               china Terrestrial Reference Frame 
# ecef_results_bds = df_bds_only.apply(calculate_ecef, axis=1, args=(3.986004418e14,7.2921150e-5))
# df_bds_only = pd.concat([df_bds_only, ecef_results_bds], axis=1)

# # Remove all rows where X, Y, or Z could not be calculated
# df_clean_bds = df_bds_only.dropna(subset=['X', 'Y', 'Z'])
# print(f"Original rows: {len(df_bds_only)}")
# print(f"Cleaned rows: {len(df_clean_bds)}")
# print(df_clean_bds[['sv', 'time', 'X', 'Y', 'Z']])


#$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

# #1. Load your file string (Galileo )
# rinex_str_gal = open(file_path_galileo).read()

# #2. Convert to DataFrame
# df = process_rinex_to_df(rinex_str_gal)

# # 3. Apply Kepler Solver
# if not df.empty:
#     df['Ek'] = df.apply(solve_kepler, axis=1)
#     print(df[['SV_ID', 'System', 'Toe', 'Ek']].head())



# # 1. Clean up rows that are missing critical data
# df = df.dropna(subset=['time', 'A', 'e', 'M0'])

# # 2. Run the calculation and store it
# ecef_results = df.apply(calculate_galileo_ecef, axis=1)

# # 3. Check if ecef_results actually has data
# if not ecef_results.empty:
#     # Use join instead of concat to ensure indices align perfectly
#     df = df.join(ecef_results)
    
#     # 4. Now verify if 'X' exists before printing
#     if 'X' in df.columns:
#         print("--- Galileo ECEF Coordinates (Meters) ---")
#         print(df[['SV_ID', 'time', 'X', 'Y', 'Z']].head())
#     else:
#         print("Error: Columns X, Y, Z were not created. Check your function return.")
# else:
#     print("Error: ecef_results is empty. Is your input DataFrame empty?")