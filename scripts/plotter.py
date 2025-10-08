import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy import stats
from scipy.stats import norm

# A small constant to add before taking a logarithm to avoid log(0) errors.
EPSILON = 1e-10

# --- Core Comparison Plotting Functions ---

def plot_signal_comparison_overlay(original, filtered, uniform_q, mu_law_q, sampling_rate, output_path):
    """Plots the entire signal in the time domain for all four stages."""
    plt.figure(figsize=(15, 7))
    time_axis = np.arange(len(original)) / sampling_rate
    plt.plot(time_axis, original, label='Original (Normalized)', alpha=0.6)
    plt.plot(time_axis, filtered, label='Filtered (1-40 Hz)', alpha=0.9, linewidth=1.5)
    
    end_index = int(2 * sampling_rate) # Show first 2 seconds
    plt.plot(time_axis[:end_index], uniform_q[:end_index], label='Uniform Q (8-bit)', alpha=1.0, linewidth=1.0, drawstyle='steps-post')
    plt.plot(time_axis[:end_index], mu_law_q[:end_index], label='μ-law Q (8-bit)', alpha=1.0, linewidth=1.0, drawstyle='steps-post', linestyle='--')

    plt.title(f"Full Signal in Time Domain Comparison\n({len(original)} Samples at {sampling_rate} Hz)")
    plt.xlabel('Time (s)'); plt.ylabel('Amplitude'); plt.grid(True)
    plt.legend(loc='lower center', ncol=4)
    plt.tight_layout(); plt.savefig(os.path.join(output_path, "comparison_signal_full.png"))

def plot_autocorrelation_comparison_overlay(original_acf, filtered_acf, uniform_q_acf, mu_law_q_acf, sampling_rate, output_path):
    max_lags = min(len(original_acf), 1024)
    lag_axis = np.arange(max_lags) / sampling_rate
    plt.figure(figsize=(15, 7))
    plt.plot(lag_axis, original_acf[:max_lags], label='Original', alpha=0.5)
    plt.plot(lag_axis, filtered_acf[:max_lags], label='Filtered', alpha=0.8)
    plt.plot(lag_axis, uniform_q_acf[:max_lags], label='Uniform Q', alpha=0.9, linestyle='--')
    plt.plot(lag_axis, mu_law_q_acf[:max_lags], label='μ-law Q', alpha=0.9, linestyle=':')
    plt.title('Direct Comparison: Autocorrelation')
    plt.xlabel('Lag (s)'); plt.ylabel('Normalized Correlation'); plt.grid(True)
    plt.legend(loc='lower center', ncol=4)
    plt.tight_layout(); plt.savefig(os.path.join(output_path, "comparison_autocorrelation.png"))

def plot_psd_comparison_overlay(psd_data, n_fft, sampling_rate, method_name, perf_data, output_path):
    freq_axis = np.fft.rfftfreq(n_fft, d=1.0 / sampling_rate)
    fig, ax = plt.subplots(figsize=(14, 8)); max_freq = 60
    mask = freq_axis <= max_freq

    ax.plot(freq_axis[mask], psd_data['orig'][mask], label='Original', alpha=0.5)
    ax.plot(freq_axis[mask], psd_data['filt'][mask], label='Filtered', linewidth=1.5)
    ax.plot(freq_axis[mask], psd_data['u_quant'][mask], label='Uniform Q', alpha=0.8, linestyle=':')
    ax.plot(freq_axis[mask], psd_data['m_quant'][mask], label='μ-law Q', alpha=0.8, linestyle='-.')
    
    ax.axvspan(1, 40, color='gray', alpha=0.2, label='Passband (1-40 Hz)')
    ax.set_title(f'Direct Comparison - {method_name}\nFour-Stage Processing Pipeline')
    ax.set_xlabel('Frequency (Hz)'); ax.set_ylabel('Power/Frequency (dB)'); ax.grid(True, which='both', linestyle='--')
    ax.legend(loc='lower center', ncol=5)

    all_min = [np.min(psd[mask]) for psd in psd_data.values()]
    all_max = [np.max(psd[mask]) for psd in psd_data.values()]
    ax.set_ylim(min(all_min) - 5, max(all_max) + 5)

    props = dict(boxstyle='round', facecolor='wheat', alpha=0.9)
    text_lines = ["Performance Metrics:\n---------------------------------"]
    stages = [('Original', ''), ('Filtered', '_filtered'), ('Uniform Q', '_uniform_quantized'), ('μ-law Q', '_mu_law_quantized')]
    for name, suffix in stages:
        perf = perf_data.get(f"{method_name}{suffix}", {})
        time_s = float(perf.get('ExecutionTime_s', 0))
        cpu_u = float(perf.get('CPUUserTime_us', 0)) / 1e6
        cpu_s = float(perf.get('CPUSystemTime_us', 0)) / 1e6
        mem = perf.get('PeakMemory_kb', 'N/A')
        text_lines.append(f"{name+':':<12} Exec: {time_s:.4f}s | CPU(U/S): {cpu_u:.4f}/{cpu_s:.4f}s | Mem: {mem}KB")
    
    ax.text(0.02, 0.98, "\n".join(text_lines), transform=ax.transAxes, fontsize=8, va='top', ha='left', bbox=props, fontfamily='monospace')

    plt.tight_layout(); plt.savefig(os.path.join(output_path, f"comparison_{method_name.lower().replace(' ', '_')}.png"))

def plot_welch_comparison_subplots(psd_data, perf_data, sampling_rate, output_path):
    window_sizes = sorted(psd_data['orig'].keys())
    fig, axes = plt.subplots(3, 2, figsize=(18, 15), sharey=False)
    fig.suptitle(f"Direct Comparison: Welch's Method\n(Four-Stage Processing Pipeline)", fontsize=16)
    axes_flat = axes.flatten()
    
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)

    for i, window_size in enumerate(window_sizes):
        ax = axes_flat[i]
        freq_axis = np.fft.rfftfreq(window_size, d=1.0 / sampling_rate)
        mask = freq_axis <= 60

        ax.plot(freq_axis[mask], psd_data['orig'][window_size][mask], label='Original', alpha=0.5)
        ax.plot(freq_axis[mask], psd_data['filt'][window_size][mask], label='Filtered', linewidth=1.5)
        ax.plot(freq_axis[mask], psd_data['u_quant'][window_size][mask], label='Uniform Q', linestyle=':')
        ax.plot(freq_axis[mask], psd_data['m_quant'][window_size][mask], label='μ-law Q', linestyle='-.')
        
        ax.axvspan(1, 40, color='gray', alpha=0.2)
        ax.set_title(f"Welch ({window_size}-point window)")
        ax.grid(True, which='both', linestyle='--')
        
        all_min = [np.min(psd[window_size][mask]) for psd in psd_data.values()]
        all_max = [np.max(psd[window_size][mask]) for psd in psd_data.values()]
        ax.set_ylim(min(all_min) - 5, max(all_max) + 5)

        text_lines = ["Perf (Exec ms | Sys ms | Mem KB)"]
        stages = [('O', ''), ('F', '_filtered'), ('UQ', '_uniform_quantized'), ('MQ', '_mu_law_quantized')]
        for name, suffix in stages:
            perf = perf_data.get(f"Welch_{window_size}{suffix}", {})
            time_ms = float(perf.get('ExecutionTime_s', 0)) * 1000
            cpu_s_ms = float(perf.get('CPUSystemTime_us', 0)) / 1000
            mem = perf.get('PeakMemory_kb', 'N/A')
            text_lines.append(f"{name+':':<3} {time_ms:6.2f} | {cpu_s_ms:6.2f} | {mem}")
        
        ax.text(0.98, 0.02, "\n".join(text_lines), transform=ax.transAxes, fontsize=7,
                va='bottom', ha='right', bbox=props, fontfamily='monospace')

    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    for ax in axes[-1, :]: ax.set_xlabel('Frequency (Hz)')
    for ax in axes[:, 0]: ax.set_ylabel('Power/Frequency (dB)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "comparison_welch_subplots.png"))

# NEW: Comprehensive 3x3 distribution plot for all methods and stages.
def plot_full_distribution_comparison(psd_data, output_path):
    """
    Generates a 3x3 grid of histograms to compare PSD distributions.
    Rows: Periodogram, Welch, Multitaper
    Columns: Filtered, Uniform Q, Mu-law Q
    """
    methods = ['Periodogram', 'Welch (1024-point)', 'Multitaper']
    stages = ['Filtered', 'Uniform Q (8-bit)', 'μ-law Q (8-bit)']
    stage_keys = ['filt', 'u_quant', 'm_quant']

    fig, axes = plt.subplots(3, 3, figsize=(20, 18), sharex=True, sharey=True)
    fig.suptitle('Comprehensive Comparison: Distribution of PSD Values (dB)', fontsize=20)

    for i, method in enumerate(methods):
        # Get all data for this method to set common bins
        all_method_data = np.concatenate([
            psd_data[method]['filt'],
            psd_data[method]['u_quant'],
            psd_data[method]['m_quant']
        ])
        bins = np.linspace(all_method_data.min(), all_method_data.max(), 150)

        for j, stage in enumerate(stages):
            ax = axes[i, j]
            stage_key = stage_keys[j]
            data = psd_data[method][stage_key]

            ax.hist(data, bins=bins, alpha=0.7, label=stage, density=True)
            
            try:
                mu, std = norm.fit(data)
                x_fit = np.linspace(data.min(), data.max(), 100)
                ax.plot(x_fit, norm.pdf(x_fit, mu, std), color='red', linestyle='--', linewidth=2, label=f'Fit (σ={std:.2f})')
            except Exception as e:
                print(f"Could not fit Gaussian for {method} - {stage}: {e}")

            ax.legend()
            ax.grid(True, linestyle=':')

            # Set titles for rows and columns
            if j == 0:
                ax.set_ylabel(f"{method}\n\nDensity", fontsize=14)
            if i == 0:
                ax.set_title(stage, fontsize=16)
            if i == len(methods) - 1:
                ax.set_xlabel('Power (dB)')

    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(output_path, "comparison_all_distributions.png"))

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
        print(f"Warning: Could not read performance file: {e}");

    try:
        # Load data for all four stages
        signal_orig = np.loadtxt(os.path.join(data_dir, "signal.txt"))
        signal_filt = np.loadtxt(os.path.join(data_dir, "signal_filtered.txt"))
        signal_u_quant = np.loadtxt(os.path.join(data_dir, "signal_uniform_quantized.txt"))
        signal_m_quant = np.loadtxt(os.path.join(data_dir, "signal_mu_law_quantized.txt"))

        acf_orig = np.loadtxt(os.path.join(data_dir, "acf.txt")); acf_orig /= acf_orig[0]
        acf_filt = np.loadtxt(os.path.join(data_dir, "acf_filtered.txt")); acf_filt /= acf_filt[0]
        acf_u_quant = np.loadtxt(os.path.join(data_dir, "acf_uniform_quantized.txt")); acf_u_quant /= acf_u_quant[0]
        acf_m_quant = np.loadtxt(os.path.join(data_dir, "acf_mu_law_quantized.txt")); acf_m_quant /= acf_m_quant[0]

        psd_periodogram = {
            'orig': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram.txt")) + EPSILON),
            'filt': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram_filtered.txt")) + EPSILON),
            'u_quant': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram_uniform_quantized.txt")) + EPSILON),
            'm_quant': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "periodogram_mu_law_quantized.txt")) + EPSILON)
        }
        psd_multitaper = {
            'orig': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper.txt")) + EPSILON),
            'filt': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper_filtered.txt")) + EPSILON),
            'u_quant': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper_uniform_quantized.txt")) + EPSILON),
            'm_quant': 10 * np.log10(np.loadtxt(os.path.join(data_dir, "multitaper_mu_law_quantized.txt")) + EPSILON)
        }
        psd_welch = {
            'orig': {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}.txt")) + EPSILON) for w in config['welch_windows']},
            'filt': {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}_filtered.txt")) + EPSILON) for w in config['welch_windows']},
            'u_quant': {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}_uniform_quantized.txt")) + EPSILON) for w in config['welch_windows']},
            'm_quant': {w: 10 * np.log10(np.loadtxt(os.path.join(data_dir, f"welch_{w}_mu_law_quantized.txt")) + EPSILON) for w in config['welch_windows']}
        }
    except IOError as e:
        print(f"Error loading data files: {e}"); return
    except (IndexError, ValueError):
        print("Error: Data appears to be empty or invalid. Cannot normalize/process. Exiting."); return

    print("\n--- Generating four-stage comparison plots ---")
    plot_signal_comparison_overlay(signal_orig, signal_filt, signal_u_quant, signal_m_quant, config['sampling_rate'], plot_dir)
    plot_autocorrelation_comparison_overlay(acf_orig, acf_filt, acf_u_quant, acf_m_quant, config['sampling_rate'], plot_dir)
    
    plot_psd_comparison_overlay(psd_periodogram, config['signal_length'], config['sampling_rate'], "Periodogram", perf_data, plot_dir)
    plot_psd_comparison_overlay(psd_multitaper, config['signal_length'], config['sampling_rate'], "Multitaper", perf_data, plot_dir)
    
    plot_welch_comparison_subplots(psd_welch, perf_data, config['sampling_rate'], plot_dir)
    
    # Assemble data for the new comprehensive distribution plot
    welch_window_for_dist = 1024
    comprehensive_psd_data = {
        'Periodogram': {
            'filt': psd_periodogram['filt'], 'u_quant': psd_periodogram['u_quant'], 'm_quant': psd_periodogram['m_quant']
        },
        'Welch (1024-point)': {
            'filt': psd_welch['filt'][welch_window_for_dist], 'u_quant': psd_welch['u_quant'][welch_window_for_dist], 'm_quant': psd_welch['m_quant'][welch_window_for_dist]
        },
        'Multitaper': {
            'filt': psd_multitaper['filt'], 'u_quant': psd_multitaper['u_quant'], 'm_quant': psd_multitaper['m_quant']
        }
    }
    plot_full_distribution_comparison(comprehensive_psd_data, plot_dir)

    print("\nAll plots have been saved to the 'results/plots' directory.")
    print("Displaying all plots simultaneously. Close all plot windows to exit the script.")
    
    plt.show()

if __name__ == "__main__":
    main()