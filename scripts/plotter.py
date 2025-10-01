import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy import stats
from scipy.stats import norm # NEW: Import norm for Gaussian fit

# Standard EEG frequency bands
EEG_BANDS = {
    'Delta': (0.5, 4),
    'Theta': (4, 8),
    'Alpha': (8, 12),
    'Beta': (12, 30),
    'Gamma': (30, 100)
}

# The individual plotting functions (plot_signal, plot_autocorrelation, etc.)
# remain the same. The changes are in plot_psd_distributions and main().

def plot_signal(data, sampling_rate, output_path):
    time_axis = np.arange(len(data)) / sampling_rate
    plt.figure(figsize=(12, 6))
    plt.plot(time_axis, data)
    plt.title(f"EEG Signal Feature (lag1_mean_0)\n({len(data)} Samples at {sampling_rate} Hz)")
    plt.xlabel("Time (s)"); plt.ylabel("Amplitude (arbitrary units)"); plt.grid(True); plt.tight_layout()
    plt.savefig(os.path.join(output_path, "01_eeg_signal.png"))

def plot_autocorrelation(data, sampling_rate, output_path):
    max_lags_samples = min(len(data), 1024)
    lag_axis_s = np.arange(max_lags_samples) / sampling_rate
    plt.figure(figsize=(12, 6))
    plt.plot(lag_axis_s, data[:max_lags_samples])
    plt.title(f"Autocorrelation Function\n({len(data)} Samples at {sampling_rate} Hz)")
    plt.xlabel("Lag (s)"); plt.ylabel("Correlation"); plt.grid(True); plt.tight_layout()
    plt.savefig(os.path.join(output_path, "02_autocorrelation.png"))

def plot_psd(psd_data, n_fft, sampling_rate, method_name, perf_metrics, output_path, filename):
    freq_axis = np.fft.rfftfreq(n_fft, d=1.0 / sampling_rate)
    plt.figure(figsize=(12, 7))
    max_freq_display = 60
    freq_mask = freq_axis <= max_freq_display
    log_psd_for_test = np.log10(psd_data[freq_mask])
    if len(log_psd_for_test) > 3:
        shapiro_stat, shapiro_p = stats.shapiro(log_psd_for_test)
        is_gaussian = "No (p < 0.05)" if shapiro_p < 0.05 else "Yes (p >= 0.05)"
    else:
        shapiro_p, is_gaussian = (0, "N/A")
    plt.semilogy(freq_axis[freq_mask], psd_data[freq_mask])
    plt.title(f"Power Spectral Density (PSD) - {method_name}\n({n_fft}-point FFT, Sampling Rate: {sampling_rate} Hz)")
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Power/Frequency (log scale)"); plt.grid(True, which='both', linestyle='--')
    handles, labels = plt.gca().get_legend_handles_labels()
    band_labels = set(labels)
    for band, (low, high) in EEG_BANDS.items():
        if high < max_freq_display:
            label = f'{band}: {low}-{high} Hz'
            if label not in band_labels:
                plt.axvspan(low, high, facecolor='gray', alpha=0.2, label=label)
                band_labels.add(label)
    plt.legend(loc='upper right', fontsize='small')
    time_s = float(perf_metrics.get('ExecutionTime_s', 0)); cpu_u_s = int(perf_metrics.get('CPUUserTime_us', 0)) / 1e6
    cpu_s_s = int(perf_metrics.get('CPUSystemTime_us', 0)) / 1e6; mem_kb = int(perf_metrics.get('PeakMemory_kb', 0))
    perf_text = (f"Performance Metrics:\n--------------------\n"
                 f"Execution Time (s): {time_s:.6f}\nCPU Time - User (s): {cpu_u_s:.6f}\n"
                 f"CPU Time - System (s): {cpu_s_s:.6f}\nPeak Memory (KB): {mem_kb}\n\n"
                 f"Gaussianity Test (on log PSD):\n--------------------\n"
                 f"Shapiro-Wilk p-value: {shapiro_p:.2e}\n"
                 f"Normally Distributed?: {is_gaussian}")
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.8)
    plt.figtext(0.65, 0.78, perf_text, fontsize=9, verticalalignment='top', bbox=props)
    plt.tight_layout(rect=[0, 0, 0.95, 1]); plt.savefig(os.path.join(output_path, filename))

def plot_psd_comparison(psd_periodogram, psd_welch_rep, psd_multitaper, signal_length, welch_nperseg_rep, sampling_rate, output_path):
    freq_axis_full = np.fft.rfftfreq(signal_length, d=1.0 / sampling_rate)
    freq_axis_welch = np.fft.rfftfreq(welch_nperseg_rep, d=1.0 / sampling_rate)
    plt.figure(figsize=(12, 7)); max_freq_display = 60
    mask_full = freq_axis_full <= max_freq_display; mask_welch = freq_axis_welch <= max_freq_display
    plt.semilogy(freq_axis_full[mask_full], psd_periodogram[mask_full], label='Periodogram', alpha=0.7, linewidth=1)
    plt.semilogy(freq_axis_welch[mask_welch], psd_welch_rep[mask_welch], label=f'Welch ({welch_nperseg_rep}-point window)', linewidth=2)
    plt.semilogy(freq_axis_full[mask_full], psd_multitaper[mask_full], label='Multitaper', alpha=0.7, linewidth=1)
    plt.title(f"Comparison of PSD Estimation Methods\n(Sampling Rate: {sampling_rate} Hz)")
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Power/Frequency (log scale)"); plt.grid(True, which='both', linestyle='--')
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(output_path, "06_psd_comparison.png"))

def plot_welch_comparison(psd_welch_data, sampling_rate, output_path):
    plt.figure(figsize=(12, 7)); max_freq_display = 60
    for window_size, psd_data in psd_welch_data.items():
        freq_axis = np.fft.rfftfreq(window_size, d=1.0 / sampling_rate)
        mask = freq_axis <= max_freq_display
        plt.semilogy(freq_axis[mask], psd_data[mask], label=f'Welch ({window_size}-point window)', alpha=0.8)
    plt.title(f"Welch's Method PSD with Varying Window Sizes\n(Sampling Rate: {sampling_rate} Hz)")
    plt.xlabel("Frequency (Hz)"); plt.ylabel("Power/Frequency (log scale)"); plt.grid(True, which='both', linestyle='--')
    plt.legend(); plt.tight_layout(); plt.savefig(os.path.join(output_path, "07_welch_comparison.png"))

# --- UPDATED FUNCTION ---
def plot_psd_distributions(psd_periodogram, psd_welch_rep, psd_multitaper, output_path):
    """Creates histograms with a Gaussian fit to show the distribution of log-scaled PSD values."""
    fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)
    fig.suptitle('Distribution of Power Spectral Density Values (log scale)', fontsize=16)

    # Data for histograms (log-scaled)
    log_psd_p = np.log10(psd_periodogram)
    log_psd_w = np.log10(psd_welch_rep)
    log_psd_m = np.log10(psd_multitaper)
    
    datasets = {
        'Periodogram': (log_psd_p, 'C0', axes[0]),
        'Welch (1024-point)': (log_psd_w, 'C1', axes[1]),
        'Multitaper': (log_psd_m, 'C2', axes[2])
    }

    for name, (data, color, ax) in datasets.items():
        # Plot the histogram
        counts, bins, _ = ax.hist(data, bins=100, alpha=0.75, color=color, density=False)
        
        # Calculate Gaussian fit
        mu, std = norm.fit(data)
        
        # Create x-axis for the fit curve
        x_fit = np.linspace(bins[0], bins[-1], 100)
        
        # Calculate the PDF and scale it to match the histogram's counts
        pdf = norm.pdf(x_fit, mu, std)
        bin_width = bins[1] - bins[0]
        scaled_pdf = pdf * len(data) * bin_width
        
        # Plot the fit
        ax.plot(x_fit, scaled_pdf, 'r--', linewidth=2, label=f'Gaussian Fit\n(μ={mu:.2f}, σ={std:.2f})')
        
        ax.set_title(f'{name} PSD Distribution')
        ax.set_xlabel('log10(Power)')
        ax.legend()
        ax.grid(True, linestyle=':')

    axes[0].set_ylabel('Count')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "08_psd_distributions.png"))

# --- UPDATED FUNCTION ---
def main():
    data_dir = os.path.join("results", "data")
    plot_dir = os.path.join("results", "plots")
    perf_file = os.path.join("results", "performance", "performance.txt")
    config_file = os.path.join(data_dir, "config.txt")
    
    print("Python plotter started...")

    # Robust config file parsing
    config = {}
    try:
        with open(config_file, 'r') as f:
            for line in f:
                parts = line.strip().split(',')
                key = parts[0]
                values = parts[1:]
                if key == 'sampling_rate': config[key] = float(values[0])
                elif key == 'signal_length': config[key] = int(values[0])
                elif key == 'welch_windows': config[key] = [int(v) for v in values]
        sampling_rate = config['sampling_rate']
        signal_length = config['signal_length']
        welch_windows = config['welch_windows']
    except (IOError, KeyError, IndexError, ValueError) as e:
        print(f"Error reading config file '{config_file}': {e}. Exiting."); return

    # Load data
    try:
        signal = np.loadtxt(os.path.join(data_dir, "signal.txt"))
        acf = np.loadtxt(os.path.join(data_dir, "acf.txt"))
        psd_periodogram = np.loadtxt(os.path.join(data_dir, "periodogram.txt"))
        psd_multitaper = np.loadtxt(os.path.join(data_dir, "multitaper.txt"))
        psd_welch_data = {w: np.loadtxt(os.path.join(data_dir, f"welch_{w}.txt")) for w in welch_windows}
    except IOError as e:
        print(f"Error loading data files from '{data_dir}': {e}"); return

    # Load performance metrics
    perf_data = {}
    try:
        with open(perf_file, 'r', newline='') as f:
            reader = csv.DictReader(f)
            for row in reader:
                perf_data[row['Method']] = row
    except (IOError, KeyError) as e:
        print(f"Error reading performance file '{perf_file}': {e}"); return

    # --- Create all plot figures ---
    print("Generating all plot figures...")
    plot_signal(signal, sampling_rate, plot_dir)
    plot_autocorrelation(acf, sampling_rate, plot_dir)
    
    if "Periodogram" in perf_data:
        print("Generating Periodogram plot...")
        plot_psd(psd_periodogram, signal_length, sampling_rate, "Periodogram", perf_data["Periodogram"], plot_dir, "03_psd_periodogram.png")

    # Loop to create an individual plot for each Welch window
    for window_size in welch_windows:
        perf_key = f"Welch_{window_size}"
        if perf_key in perf_data:
            print(f"Generating individual plot for Welch {window_size}...")
            method_title = f"Welch ({window_size}-point window)"
            file_name = f"04_psd_welch_{window_size}.png"
            plot_psd(psd_welch_data[window_size], window_size, sampling_rate, method_title, perf_data[perf_key], plot_dir, file_name)

    if "Multitaper" in perf_data:
        print("Generating Multitaper plot...")
        plot_psd(psd_multitaper, signal_length, sampling_rate, "Multitaper", perf_data["Multitaper"], plot_dir, "05_psd_multitaper.png")

    # Use a representative Welch window (e.g., 1024) for the summary plots
    rep_welch_window = 1024
    rep_welch_key = f"Welch_{rep_welch_window}"
    if rep_welch_key in perf_data:
        print("Generating general PSD comparison plot...")
        plot_psd_comparison(psd_periodogram, psd_welch_data[rep_welch_window], psd_multitaper, signal_length, rep_welch_window, sampling_rate, plot_dir)
        
        print("Generating PSD distribution plot with Gaussian fit...")
        plot_psd_distributions(psd_periodogram, psd_welch_data[rep_welch_window], psd_multitaper, plot_dir)

    print("Generating Welch window comparison plot...")
    plot_welch_comparison(psd_welch_data, sampling_rate, plot_dir)

    # --- Show all created figures at once ---
    print("Displaying all plots simultaneously. Close all windows to exit the script.")
    plt.show()

    print("All plots have been saved.")

if __name__ == "__main__":
    main()