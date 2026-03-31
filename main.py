# Section 1 — Generate a Reference Signal
#
# Choose a short DNA sequence to work with (e.g., 50–200 bases)
# Map each base or k-mer to a mean ionic current level (you can use published ONT squiggle values or define your own lookup table)
# Use NumPy to generate a clean "ground truth" current signal — a flat level per base with slight dwell-time variation
from operator import index

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sn

## Creating reference data table
df = pd.read_csv("template_median68pA.model",sep="\t",index_col=0)
template_dict = df.to_dict("index")

data = pd.DataFrame(template_dict).transpose().reset_index()
data.rename(columns={'index':'kmer'}, inplace=True)
print("Reference Data Table: \n", data)

#
dna = "ATGCGTACGTTAGCCTAGCATGCAGTCGATCGTACGATCGATGCATGCTAGCATGCGTACGATCGATCGTAGCATGCTAGCATCGATGCTAGCATGCATCGATCGTAGCATGCTAGCATCGATGCTAGCATGCATCGATCG"
#

## Finding all kmers within dna
k=6
kmers = [dna[i:i+k] for i in range(0,len(dna) - k + 1)]
kmer_df = pd.DataFrame(kmers)
kmer_df.rename(columns={0:'kmer'}, inplace=True)
print("kmer reference table: \n", kmer_df)

## creating new merged dataframe merging kmers to corresponding level_means
rg = data.merge(kmer_df,on="kmer",how='right')

## creating 1D array with dwell time 8 to mimic real nanopore levels
clean_signal = np.array(rg['level_mean'])
clean_signal = np.repeat(clean_signal,8)

# Section 2 — Add Noise
#
# Add Gaussian (thermal) noise using NumPy's random module at a defined standard deviation
# Experiment with multiple noise levels (low, medium, high) to create a set of noisy signals

np.random.seed(42)
sigma_low = 2
sigma_medium = 5
sigma_high = 10
noise_signal_low = clean_signal + np.random.normal(0,sigma_low,len(clean_signal))
noise_signal_medium = clean_signal + np.random.normal(0,sigma_medium,len(clean_signal))
noise_signal_high = clean_signal + np.random.normal(0,sigma_high,len(clean_signal))

# print("clean signal: \n", clean_signal[:25])
# print("low: \n", noise_signal_low[:25])
# print("medium: \n", noise_signal_medium[:25])
# print("high: \n", noise_signal_high[:25])

# Section 3 — Compute SNR
#
# Calculate SNR for each noise level across the signal
# Compute both global SNR (whole signal) and per-base SNR (localized)
# Store all results in a Pandas DataFrame for easy manipulation and export

def compute_snr(signal,noise_f):
    signal_power = np.mean(signal ** 2)
    noise = noise_f - signal
    noise_power = np.mean(noise ** 2)
    snr = 10 * np.log10(signal_power / noise_power)
    return snr

snr_low = compute_snr(clean_signal,noise_signal_low)
snr_medium = compute_snr(clean_signal,noise_signal_medium)
snr_high = compute_snr(clean_signal,noise_signal_high)
print("snr_low: \n", snr_low)
print("snr_medium: \n", snr_medium)
print("snr_high: \n", snr_high)

# Section 4 — Visualize with Matplotlib
#
# Plot 1: Clean signal vs. noisy signal overlay (side by side or overlaid)
# Plot 2: SNR vs. noise level (line or bar chart)
# Plot 3: Per-base SNR heatmap along the sequence

## Clean Signal vs. Noisy Signal Plot

f1 = plt.figure(1,figsize = (12,4))
plt.plot(clean_signal,label="Clean Signal",
         color = 'cyan',
         linewidth = 1)
plt.plot(noise_signal_low,
         label="Noisy Signal (Low)",
         color='gray', alpha = 0.6,
         linewidth = 0.8)
plt.xlabel("Sample Index")
plt.ylabel("Current (pA)")
plt.title("Clean vs. Noisy Signal")
plt.legend()
plt.tight_layout()
plt.savefig("Clean_vs_Noisy")

## SNR vs. Noise Level Plot
snr_values = [snr_low,snr_medium,snr_high]
labels = ["Low(σ=2)", "Medium (σ=5)", "High (σ=10)"]
f2 = plt.figure(2,figsize = (7,5))
plt.bar(labels, snr_values, color=["darkcyan",
                                   "darkturquoise",
                                   "aqua"])

plt.xlabel("Noise Level")
plt.ylabel("Signal-to-Noise Ratio")
plt.title("SNR vs. Noise Level")
plt.tight_layout()
plt.savefig("SNR_vs_Noise")


## Per-kmer SNR heatmap
dwell = 8
kmer_value = [clean_signal[i:i+dwell] for i in range(0,len(clean_signal),dwell)]
noise_value_low = [noise_signal_low[i:i+dwell] for i in range(0, len(noise_signal_low),dwell)]
noise_value_medium = [noise_signal_medium[i:i+dwell] for i in range(0, len(noise_signal_medium),dwell)]
noise_value_high = [noise_signal_high[i:i+dwell] for i in range(0, len(noise_signal_high),dwell)]

kmer_value = np.mean(np.array(kmer_value),axis=1)
noise_value_low = np.mean(np.array(noise_value_low),axis=1)
noise_value_medium = np.mean(np.array(noise_value_medium),axis=1)
noise_value_high = np.mean(np.array(noise_value_high),axis=1)

def compute_individual_snr(signal,noise_f):
    signal_power = signal ** 2
    noise = noise_f - signal
    noise_power = noise ** 2
    snr = 10 * np.log10(signal_power / noise_power)
    return snr

snr_calc = compute_individual_snr(kmer_value,noise_value_low)
# Note that this heatmap uses only low noise data. You can swap this out for medium or high noise.

snr_plt = pd.DataFrame(snr_calc).transpose().rename(columns=kmer_df["kmer"])
snr_plt.to_csv("SNR by kmer")
# print("SNR by kmer: \n", snr_plt.to_string())

f3 = plt.figure(3,figsize = (15,5))
sn.heatmap(snr_plt, cmap="RdYlGn", xticklabels=True)
plt.xticks(fontsize=4, rotation=90)
plt.xlabel("k-mer")
plt.ylabel("Signal-to-Noise (Low) Ratio")
plt.title("SNR of k-mers")
plt.tight_layout()
plt.savefig("SNR_kmers_low")

# Section 5 — Parameter Sweep & Analysis
#
# Loop over a range of noise standard deviations
# Record how SNR degrades as noise increases
# Identify a "threshold" SNR below which bases would be hard to distinguish
# Use Pandas to summarize and export results as a CSV

snr_list = np.array([])
for s in range(1,21):
    np.random.seed(12)
    noise_s = clean_signal + np.random.normal(0,s,len(clean_signal))
    noise_pd_i = [noise_s[i:i + 8] for i in range(0, len(noise_s), 8)]
    noise_pd_m = np.mean(np.array(noise_pd_i),axis=1)
    noise_pd_c = compute_snr(kmer_value,noise_pd_m)
    snr_i = np.mean(np.array(noise_pd_c))
    snr_list = np.append(snr_list,snr_i)

noise_pd_f = pd.DataFrame(snr_list).reset_index().rename(columns={0:"SNR (dB)", 'index':"sigma"})
noise_pd_f["sigma"] = noise_pd_f["sigma"] + 1
noise_pd_f.set_index("sigma", inplace=True)

noise_pd_f.to_csv("SNR by Varying Standard Deviations")
# print("SNR by varying std: \n", noise_pd_f)

# Plotting
sigmas = list(range(1,21))
f4 = plt.figure(4,figsize = (8,5))
plt.plot(noise_pd_f, "o", color="darkturquoise")
plt.xticks(sigmas)
plt.xlabel("Standard Deviation")
plt.ylabel("SNR (dB)")
plt.title("Signal-to-Noise Ratio at Varying Standard Deviations")
plt.tight_layout()
plt.savefig("SNR_20Std")

plt.show()