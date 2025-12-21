from dependencies_for_ex3 import * 



def get_precise_orbit_30s(file_path):
    """
    Reads an SP3 file, parses it, and interpolates to 30s.
    Adaptive interpolation order based on available data points.
    """
    
    # --- 1. Read File ---
    with open(file_path, 'r') as f:
            content = f.read()

    # --- 2. Parse SP3 ---
    records = []
    current_time = None
    lines = content.strip().split('\n')
    
    for line in lines:
        if line.startswith('*'):
            parts = line.split()
            # HH MM SS -> Seconds of Day
            h, m, s = int(parts[4]), int(parts[5]), float(parts[6])
            current_time = h * 3600 + m * 60 + s
        elif line.startswith('P'):
            sv_id = line[1:4].strip()
            try:
                records.append({
                    'sv': sv_id, 't_sec': current_time,
                    'X': float(line[4:18]) * 1000,
                    'Y': float(line[18:32]) * 1000,
                    'Z': float(line[32:46]) * 1000
                })
            except ValueError: continue
                
    df_raw = pd.DataFrame(records)
    if df_raw.empty:
        print("No data points found.")
        return df_raw

    # --- 3. Adaptive Interpolation ---
    target_interval = 30
    interpolated_list = []
    
    for sv_id, group in df_raw.groupby('sv'):
        group = group.sort_values('t_sec')
        n_points = len(group)
        
        # Determine the safest order (Max 9, Min 1)
        current_order = min(9, n_points - 1)
        
        if current_order < 1:
            print(f"Skipping {sv_id}: Only 1 data point available.")
            continue

        t_orig = group['t_sec'].values
        x_orig, y_orig, z_orig = group['X'].values, group['Y'].values, group['Z'].values
        
        # Grid from first to last available time
        t_target = np.arange(t_orig.min(), t_orig.max() + target_interval, target_interval)
        
        for t_t in t_target:
            # For testing with few points, use all available points
            # For real files, use a sliding window of order+1
            idx = np.abs(t_orig - t_t).argmin()
            points_needed = current_order + 1
            start = max(0, min(idx - (points_needed // 2), n_points - points_needed))
            end = start + points_needed
            
            t_win = t_orig[start:end]
            x_v = lagrange(t_win, x_orig[start:end])(t_t)
            y_v = lagrange(t_win, y_orig[start:end])(t_t)
            z_v = lagrange(t_win, z_orig[start:end])(t_t)
            
            interpolated_list.append({
                'sv': sv_id, 'time_sec': t_t, 'X': x_v, 'Y': y_v, 'Z': z_v
            })
            
    return pd.DataFrame(interpolated_list)

# --- Test ---
df_results = get_precise_orbit_30s('test_multy.SP3')
print(df_results[df_results['sv'] == 'G02'])

