from funcs_and_libs import *
#+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

# load each data and clean outlier observation  
bahrein = load_and_clean("estimation_theory/proj1/Datas/bahrein.xlsx")
chabahar = load_and_clean("estimation_theory/proj1/Datas/chabahar.xlsx")

# ectracting each dataset years and sort them by mean of each year data 
years_bahrein, means_bahrein = yearly_means(bahrein)
years_chabahar, means_chabahar = yearly_means(chabahar)


#+-+-+-+-+-+-+-+-+-+- interpolate water level for the years with NAN values+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-  

# bahrein
bahrein = expand_to_full_year_month(bahrein,'date')
bahrein_full = bahrein.sort_values(["year","month"]).reset_index(drop=True)
bahrein_full["water_level"] = bahrein_full["water_level"].interpolate(method="linear")


# chabaher
chabahar = expand_to_full_year_month(chabahar,'date')
chabahar_full = chabahar.sort_values(["year","month"]).reset_index(drop=True)
chabahar_full["water_level"] = chabahar_full["water_level"].interpolate(method="linear")


# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-
# +-+-+-+-+-+-+- part آ & part ب +-+-+-+-+-+-+-

# extracting water lavel values 
water_bahrein = bahrein_full["water_level"].values
water_chabahar = chabahar_full["water_level"].values

# Compute FFT and plot
freqs_b, mag_b , power_bahrein  = compute_fft_and_energy_spectrum(water_bahrein, sampling_rate=1)
freqs_c, mag_c , power_chabahar = compute_fft_and_energy_spectrum(water_chabahar, sampling_rate=1)


# avoid 0 powers that make log = infinity
power_bahrein = np.where(power_bahrein== 0, 1e-12, power_bahrein)
power_chabahar = np.where(power_chabahar == 0, 1e-12, power_chabahar)

#+-+-+-+-+-+- extracting the 3 frequencies with most energy for each station +-+-+-+-+-+-+-
bahrein_top3_frequencies ,bahrein_top3_power = get_most_energy_freqs(power_bahrein,freqs_b,50) 
chabahar_top3_frequencies ,chabahar_top3_power = get_most_energy_freqs(power_chabahar,freqs_c,50) 


data1 = {
    'frequency': bahrein_top3_frequencies,
    'power': bahrein_top3_power
}
top3_table_for_bahrein = pd.DataFrame(data1)

data2 = {
    'frequency': chabahar_top3_frequencies,
    'power': chabahar_top3_power
}
top3_table_for_chabahar = pd.DataFrame(data2)

for i in range(50):
    top3_table_for_bahrein.loc[i,'period(1/f)'] = 1/top3_table_for_bahrein.loc[i,'frequency'] 
    top3_table_for_chabahar.loc[i,'period(1/f)'] = 1/top3_table_for_chabahar.loc[i,'frequency'] 

print('top_table_for_bahrein:\n', top3_table_for_bahrein)
print('top_table_for_chabahar:\n', top3_table_for_chabahar)



# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-



# +-+-+-+-+-+- part - پ - NOISE ANALYSIS AND REMOVE IT +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-

reconstructed_b, fft_vals_b, fft_filtered_b = fft_full_denoise_reconstruct(water_bahrein, "Bahrein")
reconstructed_c, fft_vals_c, fft_filtered_c = fft_full_denoise_reconstruct(water_chabahar, "Chabahar")
# RMS 
rms_error_bahrein = compute_rms_error(water_bahrein, reconstructed_b)
print("RMS Error (Bahrein):", rms_error_bahrein)

rms_error_chabahar = compute_rms_error(water_chabahar, reconstructed_c)
print("RMS Error (Chabahar):", rms_error_chabahar)


freqs_bah, mag_bah, power_bah = compute_fft_and_energy_spectrum(water_bahrein, sampling_rate=1)
long_E_b  = power_bah[freqs_bah < 0.1].sum()
short_E_b = power_bah[freqs_bah >= 0.1].sum()
ratio_bahrein = long_E_b / short_E_b

freqs_cha, mag_cha, power_cha = compute_fft_and_energy_spectrum(water_chabahar, sampling_rate=1)
long_E_c  = power_cha[freqs_cha < 0.1].sum()
short_E_c = power_cha[freqs_cha >= 0.1].sum()
ratio_chabahar = long_E_c / short_E_c

print(ratio_chabahar)
print(ratio_bahrein)




top10_bahrein, spectrum_bahrein = compute_top_frequencies(water_bahrein, "Bahrein")
top10_chabahar, spectrum_chabahar = compute_top_frequencies(water_chabahar, "Chabahar")








'''
#************ plot both stations anual (meamn-year) diagram ************************************************************

#Plot Bahrein
plot_station(years_bahrein,means_bahrein,'Bahrein',marker='o',color='red')

# Plot Chabahar
plot_station(years_chabahar,means_chabahar,'Chabahar',marker='s',color='green')

#************ plot both station freq-mags diagram ****************************************************************


plot_spectrum(freqs_b,mag_b,title='frequency_spectrum bahrain')
plot_spectrum(freqs_c,mag_c,title='frequency_spectrum chabahar',color='red')



#********* plot both station log-log diagram ****************************************************************

plt.figure(figsize=(10,5))
plt.loglog(freqs_b, power_bahrein, label="Bahrein")
plt.loglog(freqs_c, power_chabahar, label="Chabahar")
plt.xlabel("Frequency (cycles per month)")
plt.ylabel("Power |F(f)|^2")
plt.title("Energy Spectrum Comparison (Log–Log)")
plt.grid(True, which="both")
plt.legend()
plt.show(block=False)  # <-- non-blocking



'''
# +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-






























# Keep the plots open
plt.pause(0.1)   # refresh GUI
input("Press Enter to close all figures...")  # wait until user is ready
plt.close('all')







