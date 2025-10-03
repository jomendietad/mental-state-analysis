import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy import stats
from scipy.stats import norm

# A small constant to add before taking a logarithm to avoid log(0) errors.
EPSILON = 1e-10

# --- Core Comparison Plotting Functions ---

def plot_signal_comparison_overlay(original, filtered, sampling_rate, output_path):
    """Plots the entire signal in the time domain."""
    plt.figure(figsize=(15, 7))
    time_axis = np.arange(len(original)) / sampling_rate
    plt.plot(time_axis, original, label='Original', alpha=0.7)
    plt.plot(time_axis, filtered, label='Filtered (1-40 Hz)', alpha=0.9, linewidth=1.5)
    plt.title(f"Full Signal in Time Domain\n({len(original)} Samples at {sampling_rate} Hz)")
    plt.xlabel('Time (s)'); plt.ylabel('Amplitude'); plt.grid(True)
    plt.legend(loc='lower center', ncol=2)
    plt.tight_layout(); plt.savefig(os.path.join(output_path, "comparison_signal_full.png"))
    # We don't close the plot here so plt.show() can display it

def plot_autocorrelation_comparison_overlay(original_acf, filtered_acf, sampling_rate, output_path):
    max_lags = min(len(original_acf), 1024)
    lag_axis = np.arange(max_lags) / sampling_rate
    plt.figure(figsize=(15, 7))
    plt.plot(lag_axis, original_acf[:max_lags], label='Original', alpha=0.7)
    plt.plot(lag_axis, filtered_acf[:max_lags], label='Filtered (1-40 Hz)', alpha=0.9, linewidth=1.5)
    plt.title('Direct Comparison: Autocorrelation')
    plt.xlabel('Lag (s)'); plt.ylabel('Normalized Correlation'); plt.grid(True)
    plt.legend(loc='lower center', ncol=2)
    plt.tight_layout(); plt.savefig(os.path.join(output_path, "comparison_autocorrelation.png"))
    # We don't close the plot here

def plot_psd_comparison_overlay(original_psd, filtered_psd, n_fft, sampling_rate, method_name, perf_orig, perf_filt, output_path):
    freq_axis = np.fft.rfftfreq(n_fft, d=1.0 / sampling_rate)
    fig, ax = plt.subplots(figsize=(14, 8)); max_freq = 60
    
    mask = freq_axis <= max_freq

    ax.plot(freq_axis[mask], original_psd[mask], label='Original', alpha=0.85)
    ax.plot(freq_axis[mask], filtered_psd[mask], label='Filtered', linewidth=1.5)
    
    ax.axvspan(1, 40, color='gray', alpha=0.2, label='Passband (1-40 Hz)')
    ax.set_title(f'Direct Comparison - {method_name}\nOriginal vs. Filtered (1-40 Hz)\n(Sampling Rate: {sampling_rate} Hz, FFT Points: {n_fft})', fontsize=14)
    ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel('Power/Frequency (dB)'); ax.grid(True, which='both', linestyle='--')
    ax.legend(loc='lower center', ncol=3)

    min_val = min(np.min(original_psd[mask]), np.min(filtered_psd[mask]))
    max_val = max(np.max(original_psd[mask]), np.max(filtered_psd[mask]))
    ax.set_ylim(min_val - 5, max_val + 5)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
    time_s_orig = float(perf_orig.get('ExecutionTime_s', 0)); cpu_u_s_orig = int(perf_orig.get('CPUUserTime_us', 0)) / 1e6
    text_orig = (f"Original Data:\n--------------------\n"
                 f"Exec Time: {time_s_orig:.5f} s\nCPU Time: {cpu_u_s_orig:.5f} s\n"
                 f"Peak Mem: {perf_orig.get('PeakMemory_kb', 'N/A')} KB")
    ax.text(0.02, 0.98, text_orig, transform=ax.transAxes, fontsize=8, va='top', ha='left', bbox=props)
    time_s_filt = float(perf_filt.get('ExecutionTime_s', 0)); cpu_u_s_filt = int(perf_filt.get('CPUUserTime_us', 0)) / 1e6
    text_filt = (f"Filtered Data:\n--------------------\n"
                 f"Exec Time: {time_s_filt:.5f} s\nCPU Time: {cpu_u_s_filt:.5f} s\n"
                 f"Peak Mem: {perf_filt.get('PeakMemory_kb', 'N/A')} KB")
    ax.text(0.98, 0.98, text_filt, transform=ax.transAxes, fontsize=8, va='top', ha='right', bbox=props)

    plt.tight_layout(); plt.savefig(os.path.join(output_path, f"comparison_{method_name.lower().replace(' ', '_')}.png"))
    # We don't close the plot here

def plot_welch_comparison_subplots(psd_welch_orig, psd_welch_filt, perf_data, sampling_rate, output_path):
    window_sizes = sorted(psd_welch_orig.keys())
    fig, axes = plt.subplots(3, 2, figsize=(18, 15), sharey=False)
    fig.suptitle(f"Direct Comparison: Welch's Method\n(Sampling Rate: {sampling_rate} Hz)", fontsize=16)
    axes_flat = axes.flatten()
    for i, window_size in enumerate(window_sizes):
        ax = axes_flat[i]; freq_axis = np.fft.rfftfreq(window_size, d=1.0 / sampling_rate)
        
        mask = freq_axis <= 60

        ax.plot(freq_axis[mask], psd_welch_orig[window_size][mask], label='Original', alpha=0.7)
        ax.plot(freq_axis[mask], psd_welch_filt[window_size][mask], label='Filtered', linewidth=1.5)
        
        ax.axvspan(1, 40, color='gray', alpha=0.2)
        ax.set_title(f"Welch ({window_size}-point window)"); ax.grid(True, which='both', linestyle='--')
        
        min_val = min(np.min(psd_welch_orig[window_size][mask]), np.min(psd_welch_filt[window_size][mask]))
        max_val = max(np.max(psd_welch_orig[window_size][mask]), np.max(psd_welch_filt[window_size][mask]))
        ax.set_ylim(min_val - (0.1 * abs(min_val)), max_val + (0.1 * abs(max_val)))

    handles, labels = axes_flat[0].get_legend_handles_labels(); fig.legend(handles, labels, loc='upper right')
    for ax in axes[-1, :]: ax.set_xlabel('Frequency (Hz)')
    for ax in axes[:, 0]: ax.set_ylabel('Power/Frequency (dB)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95]); plt.savefig(os.path.join(output_path, "comparison_welch_subplots.png"))
    # We don't close the plot here

def plot_distribution_comparison(psd_periodogram_orig, psd_periodogram_filt, psd_welch_rep_orig, psd_welch_rep_filt, psd_multitaper_orig, psd_multitaper_filt, output_path):
    fig, axes = plt.subplots(1, 3, figsize=(20, 7), sharey=True)
    fig.suptitle('Direct Comparison: Distribution of PSD Values (dB)', fontsize=16)
    datasets = {'Periodogram': (psd_periodogram_orig, psd_periodogram_filt, axes[0]), 'Welch (1024-point)': (psd_welch_rep_orig, psd_welch_rep_filt, axes[1]), 'Multitaper': (psd_multitaper_orig, psd_multitaper_filt, axes[2])}
    for name, (data_orig, data_filt, ax) in datasets.items():
        ax.hist(data_orig, bins=100, alpha=0.7, label='Original', density=True)
        ax.hist(data_filt, bins=100, alpha=0.7, label='Filtered', density=True)
        mu_orig, std_orig = norm.fit(data_orig); x_fit_orig = np.linspace(data_orig.min(), data_orig.max(), 100)
        ax.plot(x_fit_orig, norm.pdf(x_fit_orig, mu_orig, std_orig), color='blue', linestyle='--', linewidth=2, label=f'Original Fit (σ={std_orig:.2f})')
        mu_filt, std_filt = norm.fit(data_filt); x_fit_filt = np.linspace(data_filt.min(), data_filt.max(), 100)
        ax.plot(x_fit_filt, norm.pdf(x_fit_filt, mu_filt, std_filt), color='red', linestyle='--', linewidth=2, label=f'Filtered Fit (σ={std_filt:.2f})')
        ax.set_title(name); ax.set_xlabel('Power (dB)'); ax.legend(); ax.grid(True, linestyle=':')
    axes[0].set_ylabel('Density'); plt.tight_layout(rect=[0, 0.03, 1, 0.95]); plt.savefig(os.path.join(output_path, "comparison_distributions.png"))
    # We don't close the plot here

def main():
    data_dir = os.path.join("results", "data"); plot_dir = os.path.join("results", "plots")
    print("Python plotter started...")
    try:
        with open(os.path.join(data_dir, "config.txt"), 'r') as f:
            config_reader = csv.reader(f); config = {row[0]: [v for v in row[1:]] for row in config_reader}
        config['sampling_rate'] = float(config['sampling_rate'][0]); config['signal_length'] = int(config['signal_length'][0])
        config['welch_windows'] = [int(v) for v in config['welch_windows']]
    except Exception as e:
        print(f"Error reading config file: {e}. Exiting."); return
    perf_data = {}
    try:
        with open(os.path.join("results", "performance", "performance.txt"), 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader: perf_data[row['Method']] = row
    except Exception as e:
        print(f"Error reading performance file: {e}"); return
    try:
        signal_orig = np.loadtxt(os.path.join(data_dir, "signal.txt")); signal_filt = np.loadtxt(os.path.join(data_dir, "signal_filtered.txt"))
        acf_orig = np.loadtxt(os.path.join(data_dir, "acf.txt")); acf_orig = acf_orig / acf_orig[0]
        acf_filt = np.loadtxt(os.path.join(data_dir, "acf_filtered.txt")); acf_filt = acf_filt / acf_filt[0]
        periodogram_orig = 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram.txt")) + EPSILON)
        periodogram_filt = 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram_filtered.txt")) + EPSILON)
        multitaper_orig = 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper.txt")) + EPSILON)
        multitaper_filt = 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper_filtered.txt")) + EPSILON)
        psd_welch_orig = {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}.txt")) + EPSILON) for w in config['welch_windows']}
        psd_welch_filt = {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}_filtered.txt")) + EPSILON) for w in config['welch_windows']}
    except IOError as e:
        print(f"Error loading data files: {e}"); return
    except IndexError:
        print("Error: Autocorrelation data appears to be empty. Cannot normalize. Exiting."); return

    print("\n--- Generating direct comparison plots (Original vs. Filtered) ---")
    plot_signal_comparison_overlay(signal_orig, signal_filt, config['sampling_rate'], plot_dir)
    # MODIFIED: The call to the zoom plot function has been removed.
    plot_autocorrelation_comparison_overlay(acf_orig, acf_filt, config['sampling_rate'], plot_dir)
    plot_psd_comparison_overlay(periodogram_orig, periodogram_filt, config['signal_length'], config['sampling_rate'], "Periodogram", perf_data.get("Periodogram", {}), perf_data.get("Periodogram_filtered", {}), plot_dir)
    plot_psd_comparison_overlay(multitaper_orig, multitaper_filt, config['signal_length'], config['sampling_rate'], "Multitaper", perf_data.get("Multitaper", {}), perf_data.get("Multitaper_filtered", {}), plot_dir)
    plot_welch_comparison_subplots(psd_welch_orig, psd_welch_filt, perf_data, config['sampling_rate'], plot_dir)
    plot_distribution_comparison(periodogram_orig, periodogram_filt, psd_welch_orig[1024], psd_welch_filt[1024], multitaper_orig, multitaper_filt, plot_dir)

    print("\nAll plots have been saved to the 'results/plots' directory.")
    print("Displaying all plots simultaneously. Close all plot windows to exit the script.")
    
    plt.show()

if __name__ == "__main__":
    main()