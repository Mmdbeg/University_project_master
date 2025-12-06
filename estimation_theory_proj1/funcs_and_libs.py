# essential libraries +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

import pandas as pd 
import matplotlib.pyplot as plt
import numpy as np
from datetime import datetime, timedelta
import calendar

# essential functions +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
def outlayer_remover(df, data_col_name):
    a =df[data_col_name]
    Q1 = a.quantile(0.25)
    Q3 = a.quantile(0.75)
    IQR = Q3 - Q1

    lower = Q1 - 1.5 * IQR
    upper = Q3 + 1.5 * IQR
    # Keep only non-outliers
    cleaned_df = df[(a >= lower) & (a <= upper)]

    #print(len(cleaned_df)-len(df))
    return cleaned_df


# ----- Load and clear data from outlier -----
def load_and_clean(path):
    df = pd.read_excel(path)
    df = outlayer_remover(df, 'water_level')
    return df.reset_index(drop=True)


# -------------------- Compute Yearly Means --------------------
def yearly_means(df):
    df["year"] = df["date"].astype(float).astype(int)
    year_means = df.groupby("year")["water_level"].mean()
    return year_means.index.tolist(), year_means.tolist()


# -------------------- Plot Function --------------------
def plot_station(years, means, title ,marker='o', color='blue'):
    plt.figure(figsize=(10,5))
    plt.plot(years, means, marker=marker, label=title, color=color)
    plt.xlabel("Year")
    plt.ylabel("Mean Water Level")
    plt.title(f"Mean Yearly Water Level - {title}")
    plt.legend()
    plt.grid(True)
    plt.show(block=False)

# ----------compute FFT transformation ----------------------------

def compute_fft_and_energy_spectrum(values, sampling_rate=1):
    """
    Compute FFT magnitudes and energy spectrum |F(f)|^2.
    Returns only positive frequencies.
    """

    # Compute FFT
    fft_raw = np.fft.fft(values)
    freqs_raw = np.fft.fftfreq(len(values), d=sampling_rate)

    # Keep only positive frequencies
    mask = freqs_raw > 0
    freqs = freqs_raw[mask]
    fft_vals = np.abs(fft_raw[mask])       # magnitude spectrum
    power = fft_vals ** 2                  # energy spectrum |F(f)|^2

    return freqs, fft_vals, power


# ---------------- plot spectrum -----------------------------
def plot_spectrum(freqs,values, title="Frequency Spectrum", color='blue',):
    plt.figure(figsize=(10,5))
    plt.plot(freqs, values,color=color)
    plt.xlabel("Frequency (cycles per unit time)")
    plt.ylabel("Magnitude")
    plt.title(title)
    plt.grid(True)
    plt.show(block=False)  # <-- non-blocking


#+-+-+- create a data frame with nan for the mosing data +-+-+-+-+-+-+-+-+-+-
def expand_to_full_year_month(df, date_col="date", value_col="water_level"):
    """
    df: original dataframe with decimal-year column
    date_col: name of the decimal-year column
    value_col: name of the observation column
    Returns: dataframe with columns ['year','month',value_col] where each year has 12 months.
             Missing months have NaN in value_col.
    """
    df = df.copy()
    
    # Step 1: Convert decimal year to year and month
    df["year"] = df[date_col].astype(int)
    df["fraction"] = df[date_col] - df["year"]
    df["month"] = (df["fraction"] * 12).apply(np.floor).astype(int) + 1
    df["month"] = df["month"].clip(1, 12)
    
    # Step 2: Keep only necessary columns
    df = df[["year","month", value_col]]
    
    # Step 3: Create full year-month grid
    all_years = np.arange(df["year"].min(), df["year"].max()+1)
    all_months = np.arange(1,13)
    full_grid = pd.MultiIndex.from_product([all_years, all_months], names=["year","month"]).to_frame(index=False)
    
    # Step 4: Merge original data onto full grid
    df_full = pd.merge(full_grid, df, on=["year","month"], how="left")
    
    return df_full


# +-+-+- Sort frequencies by power (energy) and return the top k. +-+-+-+- 
def get_most_energy_freqs(power,freqs,k):
        """
    Sort frequencies by power (energy) and return the top k.

    Parameters:
        freqs : array-like
        power : array-like (P(f) = |F(f)|^2)
        k     : number of top frequencies to extract

    Returns:
        top_freqs  : array of top-k frequencies
        top_powers : their corresponding power values
    """
        freqs = np.array(freqs)
        power = np.array(power)
        # Sort indices by descending power
        sorted_idx = np.argsort(power)[::-1]

        top_freqs = freqs[sorted_idx][:k]
        top_powers = power[sorted_idx][:k]



        return top_freqs, top_powers




def fft_full_denoise_reconstruct(values, station_name, threshold=0.005):

    # 1. FFT
    fft_vals = np.fft.fft(values)

    # 2. Power spectrum
    power = np.abs(fft_vals)**2

    # ---- NEW: ignore DC (index 0) when finding max power ----
    power_no_dc = power[1:]               # drop element 0
    max_non_dc = power_no_dc.max()
    
    # threshold based on *dominant non-DC frequency*
    lower_bound = max_non_dc * threshold
    # ---------------------------------------------------------

    # 3. Filter FFT
    fft_filtered = fft_vals.copy()
    fft_filtered[power < lower_bound] = 0

    # 4. Keep DC term unchanged (optional but recommended)
    fft_filtered[0] = fft_vals[0]

    # 5. iFFT → return real signal
    reconstructed = np.fft.ifft(fft_filtered).real

    # 6. Plot
    plt.figure(figsize=(10,4))
    plt.plot(values, label="Original Signal", alpha=0.8)
    plt.plot(reconstructed, label="Denoised Signal")
    plt.legend()
    plt.title(f"Before & After FFT Denoising ({station_name})")
    plt.xlabel("Time index")
    plt.ylabel("Value")
    plt.grid(True)
    plt.show(block=False)

    return reconstructed, fft_vals, fft_filtered


def compute_rms_error(original, reconstructed):
    diff = original - reconstructed
    rms = np.sqrt(np.mean(diff**2))
    return rms


def make_denoised_freq_dataframe(values, fft_filtered):
    # Number of samples
    N = len(values)
    sampling_rate = 1  # monthly sampling

    # Frequency bins
    freqs = np.fft.fftfreq(N, d=sampling_rate)

    # Mask: keep only non-zero (denoised) frequencies
    mask = fft_filtered != 0

    # Build dataframe
    df_denoised = pd.DataFrame({
        "frequency": freqs[mask],
        "power": (np.abs(fft_filtered[mask]))**2
    })

    # Sort by power descending
    df_denoised = df_denoised.sort_values(by="power", ascending=False).reset_index(drop=True)

    return df_denoised

def make_single_side_spectrum(values, fft_filtered):
    N = len(values)
    freqs = np.fft.fftfreq(N, d=1)       # sampling=1 month
    power = np.abs(fft_filtered)**2

    df = pd.DataFrame({
        "frequency": freqs,
        "power": power
    })

    # 1) Keep only positive and zero frequencies
    df = df[df["frequency"] >= 0].copy()

    # 2) Double the energy except DC & Nyquist
    df.loc[df["frequency"] > 0, "power"] *= 2

    # If N is even → include Nyquist freq
    # If N is odd → there is no exact Nyquist freq
    return df.sort_values("frequency").reset_index(drop=True)



def compute_top_frequencies(values, station_name):
    
    N = len(values)
    dt = 1  # 1 month sampling
    fft_vals = np.fft.fft(values)
    freqs = np.fft.fftfreq(N, d=dt)

    # Only positive frequencies
    mask = freqs > 0
    freqs = freqs[mask]
    power = np.abs(fft_vals[mask])**2

    # Build dataframe
    df = pd.DataFrame({"frequency": freqs, "power": power})

    # Sort by highest power
    df_top10 = df.sort_values(by="power", ascending=False).head(10).reset_index(drop=True)

    print(f"\n==============================")
    print(f"   TOP 10 FREQUENCIES — {station_name}")
    print(f"==============================")
    print(df_top10)

    # Plot full spectrum + mark top 10
    plt.figure(figsize=(10,5))
    plt.plot(freqs, power, alpha=0.6, label="Power Spectrum")
    plt.scatter(df_top10["frequency"], df_top10["power"], color='red', label="Top 10 Peaks")
    
    for i in range(10):
        f = df_top10.loc[i, "frequency"]
        p = df_top10.loc[i, "power"]
        plt.text(f, p, f"{i+1}", fontsize=8, color='red')

    plt.xlabel("Frequency (cycles/month)")
    plt.ylabel("Power (|FFT|²)")
    plt.title(f"Frequency Spectrum — {station_name}")
    plt.grid(alpha=0.3)
    plt.legend()
    plt.show(block=False)

    return df_top10, df