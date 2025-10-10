import numpy as np
import matplotlib.pyplot as plt
import csv
import os
from scipy import signal as sp_signal
from scipy import stats
from scipy.stats import norm

EPSILON = 1e-10

# --- NEW PLOTTING FUNCTIONS AS REQUESTED ---

def plot_time_domain_comparison(original, filtered, uniform_q, mu_law_q, sampling_rate, output_path):
    fig, axes = plt.subplots(2, 1, figsize=(15, 10), sharex=True)
    fig.suptitle('New Request: Time Domain Signal Comparison', fontsize=16)
    time_axis = np.arange(len(original)) / sampling_rate
    axes[0].plot(time_axis, original, label='Original (Normalized)', alpha=0.7)
    axes[0].plot(time_axis, filtered, label='Filtered', alpha=1.0, linewidth=1.2)
    axes[0].set_title('Original vs. Filtered Signal')
    axes[0].set_ylabel('Amplitude')
    axes[0].grid(True); axes[0].legend()
    axes[1].plot(time_axis, uniform_q, label='Uniform Quantization (8-bit)', linewidth=1.0)
    axes[1].plot(time_axis, mu_law_q, label='μ-law Quantization (8-bit)', linestyle='--', linewidth=1.0)
    axes[1].set_title('Comparison of Quantized Signals')
    axes[1].set_xlabel('Time (s)'); axes[1].set_ylabel('Amplitude')
    axes[1].grid(True); axes[1].legend()
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "new_time_domain_comparison.png"))

def plot_new_acf_comparison(acf_data, sampling_rate, output_path):
    plt.figure(figsize=(15, 8))
    max_lags = min(len(acf_data['filt']), 4096)
    lag_axis = np.arange(max_lags) / sampling_rate
    acf_filt = acf_data['filt'] / acf_data['filt'][0]
    acf_u_quant = acf_data['u_quant'] / acf_data['u_quant'][0]
    acf_m_quant = acf_data['m_quant'] / acf_data['m_quant'][0]
    plt.plot(lag_axis, acf_filt[:max_lags], label='Filtered', alpha=0.9)
    plt.plot(lag_axis, acf_u_quant[:max_lags], label='Uniform Q', linestyle='--', alpha=0.8)
    plt.plot(lag_axis, acf_m_quant[:max_lags], label='μ-law Q', linestyle=':', alpha=0.8)
    plt.title('New Request: Autocorrelation Comparison (Unbiased FFT Method)')
    plt.xlabel('Lag (s)'); plt.ylabel('Normalized Autocorrelation')
    plt.grid(True, which='both', linestyle='--')
    plt.legend()
    plt.xlim(0, max_lags / sampling_rate)
    plt.tight_layout()
    plt.savefig(os.path.join(output_path, "new_acf_comparison.png"))

def estimate_psd_wavelet_py(signal_data, fs, num_scales=128):
    wavelet = 'morl'
    freqs = np.logspace(np.log10(0.5), np.log10(fs/2), num=num_scales)
    scales = (6 / (2 * np.pi)) * (fs / freqs)
    cwt_matrix = sp_signal.cwt(signal_data, sp_signal.morlet2, scales, w=6)
    power = np.abs(cwt_matrix)**2
    psd = np.mean(power, axis=1)
    return freqs, psd

def plot_new_psd_comparison(psd_data, sampling_rate, output_path):
    methods = ['Periodogram', 'Welch', 'Multitaper', 'Wavelet']
    fig, axes = plt.subplots(len(methods), 1, figsize=(12, 18), sharex=True)
    fig.suptitle('New Request: Power Spectral Density Comparison Across All Methods', fontsize=16)
    for i, method in enumerate(methods):
        ax = axes[i]
        data = psd_data[method]
        ax.semilogx(data['freqs_filt'], data['psd_filt'], label='Filtered')
        ax.semilogx(data['freqs_u_quant'], data['psd_u_quant'], label='Uniform Q', linestyle='--')
        ax.semilogx(data['freqs_m_quant'], data['psd_m_quant'], label='μ-law Q', linestyle=':')
        ax.set_title(method)
        ax.set_ylabel('Power/Frequency (dB)')
        ax.grid(True, which='both', linestyle='--')
        ax.legend()
        ax.set_xlim(0.5, sampling_rate / 2)
    axes[-1].set_xlabel('Frequency (Hz) [Log Scale]')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(output_path, "new_psd_comparison.png"))

# --- ORIGINAL PLOTTING FUNCTIONS (RESTORED AND UPDATED) ---

def plot_original_psd_comparison(psd_data, n_fft, sampling_rate, method_name, perf_data, output_path):
    freq_axis = np.fft.rfftfreq(n_fft, d=1.0 / sampling_rate)
    fig, ax = plt.subplots(figsize=(14, 8)); max_freq = 60
    mask = freq_axis <= max_freq
    ax.plot(freq_axis[mask], psd_data['orig'][mask], label='Original', alpha=0.5)
    ax.plot(freq_axis[mask], psd_data['filt'][mask], label='Filtered', linewidth=1.5)
    ax.plot(freq_axis[mask], psd_data['u_quant'][mask], label='Uniform Q (8-bit)', alpha=0.8, linestyle=':')
    ax.plot(freq_axis[mask], psd_data['m_quant'][mask], label='μ-law Q (8-bit)', alpha=0.8, linestyle='-.')
    ax.axvspan(1, 40, color='gray', alpha=0.2, label='Passband (1-40 Hz)')
    ax.set_title(f'Original Power Spectral Density Comparison - {method_name}')
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
    plt.tight_layout(); plt.savefig(os.path.join(output_path, f"original_psd_{method_name.lower()}.png"))

def plot_welch_comparison_subplots(psd_data, perf_data, sampling_rate, output_path):
    window_sizes = sorted(psd_data['orig'].keys())
    fig, axes = plt.subplots(3, 2, figsize=(18, 15), sharey=False)
    fig.suptitle(f"Original PSD Comparison: Welch's Method", fontsize=16)
    axes_flat = axes.flatten()
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.7)
    for i, window_size in enumerate(window_sizes):
        ax = axes_flat[i]
        freq_axis = np.fft.rfftfreq(window_size, d=1.0 / sampling_rate)
        mask = freq_axis <= 60
        ax.plot(freq_axis[mask], psd_data['orig'][window_size][mask], label='Original', alpha=0.5)
        ax.plot(freq_axis[mask], psd_data['filt'][window_size][mask], label='Filtered', linewidth=1.5)
        ax.plot(freq_axis[mask], psd_data['u_quant'][window_size][mask], label='Uniform Q (8-bit)', linestyle=':')
        ax.plot(freq_axis[mask], psd_data['m_quant'][window_size][mask], label='μ-law Q (8-bit)', linestyle='-.')
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
            text_lines.append(f"{name+':':<5} {time_ms:6.2f} | {cpu_s_ms:6.2f} | {mem}")
        ax.text(0.98, 0.02, "\n".join(text_lines), transform=ax.transAxes, fontsize=7,
                va='bottom', ha='right', bbox=props, fontfamily='monospace')
    handles, labels = axes_flat[0].get_legend_handles_labels()
    fig.legend(handles, labels, loc='upper right')
    for ax in axes[-1, :]: ax.set_xlabel('Frequency (Hz)')
    for ax in axes[:, 0]: ax.set_ylabel('Power/Frequency (dB)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "original_welch_subplots.png"))

# MODIFIED: Removed shared Y-axis to allow individual scaling.
def plot_full_distribution_comparison(psd_data, output_path):
    methods = ['Periodogram', 'Welch (1024-point)', 'Multitaper']
    stages = ['Filtered', 'Uniform Quantization', 'μ-law Quantization']
    stage_keys = ['filt', 'u_quant', 'm_quant']
    fig, axes = plt.subplots(3, 3, figsize=(20, 18), sharex=True) # Removed sharey=True
    fig.suptitle('Original Comprehensive Comparison: Distribution of PSD Values (dB)', fontsize=20)
    for i, method in enumerate(methods):
        all_method_data = np.concatenate([psd_data[method][key] for key in stage_keys])
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
            except Exception: pass
            ax.legend()
            ax.grid(True, linestyle=':')
            if j == 0: ax.set_ylabel(f"{method}\n\nDensity", fontsize=14)
            if i == 0: ax.set_title(stage, fontsize=16)
            if i == len(methods) - 1: ax.set_xlabel('Power (dB)')
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])
    plt.savefig(os.path.join(output_path, "original_all_distributions.png"))

def plot_pcm_pulse_train(quantized_signal, encoded_stream, sampling_rate, output_path):
    num_samples_to_show = 32
    fig, axes = plt.subplots(2, 1, figsize=(15, 8), sharex=True, gridspec_kw={'height_ratios': [2, 3]})
    fig.suptitle(f'PCM Encoding Visualization: Serial Pulse Train (First {num_samples_to_show} Samples)', fontsize=16)
    bits_per_sample = 8
    total_bits = num_samples_to_show * bits_per_sample
    sample_period = 1.0 / sampling_rate
    bit_period = sample_period / bits_per_sample
    time_axis_samples = np.arange(num_samples_to_show) * sample_period
    axes[0].plot(time_axis_samples, quantized_signal[:num_samples_to_show], 'o-', label='Quantized Sample Value', color='C2', markersize=8)
    axes[0].set_ylabel('Amplitude'); axes[0].grid(True); axes[0].legend()
    axes[0].set_title('Input: Uniformly Quantized Signal')
    bit_stream_str = "".join(encoded_stream[:num_samples_to_show])
    bit_stream_int = [int(bit) for bit in bit_stream_str]
    time_axis_bits = np.arange(total_bits + 1) * bit_period
    plot_bit_stream = np.repeat(bit_stream_int, 2)
    plot_time_axis = np.repeat(time_axis_bits, 2)[1:-1]
    axes[1].plot(plot_time_axis, plot_bit_stream, color='C0')
    axes[1].set_ylim(-0.1, 1.1)
    axes[1].set_yticks([0, 1]); axes[1].set_yticklabels(['0', '1'])
    axes[1].set_xlabel('Time (s)'); axes[1].set_ylabel('Bit Value')
    axes[1].grid(True)
    axes[1].set_title('Output: 8-bit Serial PCM Pulse Train')
    for i in range(num_samples_to_show):
        sample_start_time = i * sample_period
        axes[0].axvline(x=sample_start_time, color='gray', linestyle='--', linewidth=0.8)
        axes[1].axvline(x=sample_start_time, color='gray', linestyle='--', linewidth=0.8)
        axes[1].text(sample_start_time + sample_period / 2, 1.2, encoded_stream[i], ha='center', va='bottom', fontsize=8)
    axes[1].axvline(x=num_samples_to_show * sample_period, color='gray', linestyle='--', linewidth=0.8)
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "visualization_pcm_pulse_train.png"))

def plot_quantization_error(filtered_signal, q8, q4, q2, sampling_rate, output_path):
    fig, axes = plt.subplots(4, 1, figsize=(15, 12), sharex=True, sharey=True)
    fig.suptitle('Quantization Error vs. Number of Bits', fontsize=16)
    time_axis = np.arange(len(filtered_signal)) / sampling_rate
    axes[0].plot(time_axis, filtered_signal, label='Reference Signal', color='black')
    axes[0].set_title('Original Filtered Signal (Input to Quantizer)')
    axes[0].set_ylabel('Amplitude'); axes[0].grid(True); axes[0].legend()
    error_8bit = filtered_signal - q8
    mse_8bit = np.mean(error_8bit**2)
    axes[1].plot(time_axis, error_8bit, label='8-bit Error', color='C0')
    axes[1].set_title(f'8-bit Quantization Error (MSE = {mse_8bit:.2e})')
    axes[1].set_ylabel('Error'); axes[1].grid(True); axes[1].legend()
    error_4bit = filtered_signal - q4
    mse_4bit = np.mean(error_4bit**2)
    axes[2].plot(time_axis, error_4bit, label='4-bit Error', color='C1')
    axes[2].set_title(f'4-bit Quantization Error (MSE = {mse_4bit:.2e})')
    axes[2].set_ylabel('Error'); axes[2].grid(True); axes[2].legend()
    error_2bit = filtered_signal - q2
    mse_2bit = np.mean(error_2bit**2)
    axes[3].plot(time_axis, error_2bit, label='2-bit Error', color='C3')
    axes[3].set_title(f'2-bit Quantization Error (MSE = {mse_2bit:.2e})')
    axes[3].set_xlabel('Time (s)'); axes[3].set_ylabel('Error')
    axes[3].grid(True); axes[3].legend()
    plt.tight_layout(rect=[0, 0.03, 1, 0.95])
    plt.savefig(os.path.join(output_path, "analysis_quantization_error.png"))

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
        # Load all signals
        signal_orig = np.loadtxt(os.path.join(data_dir, "signal.txt"))
        signal_filt = np.loadtxt(os.path.join(data_dir, "signal_filtered.txt"))
        signal_u_quant = np.loadtxt(os.path.join(data_dir, "signal_uniform_quantized.txt"))
        signal_m_quant = np.loadtxt(os.path.join(data_dir, "signal_mu_law_quantized.txt"))
        signal_q4 = np.loadtxt(os.path.join(data_dir, "signal_quantized_4bit.txt"))
        signal_q2 = np.loadtxt(os.path.join(data_dir, "signal_quantized_2bit.txt"))
        pcm_stream = np.loadtxt(os.path.join(data_dir, "pcm_encoded_stream.txt"), dtype=str)

        # Load new FFT-based ACFs
        acf_data = {
            'filt': np.loadtxt(os.path.join(data_dir, "acf_fft_filtered.txt")),
            'u_quant': np.loadtxt(os.path.join(data_dir, "acf_fft_uniform_quantized.txt")),
            'm_quant': np.loadtxt(os.path.join(data_dir, "acf_fft_mu_law_quantized.txt"))
        }

        # Load all PSDs from C for all stages
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
    except (IOError, ValueError) as e:
        print(f"Error loading data files from C: {e}"); return

    # --- Prepare data for new PSD comparison plot ---
    psd_comparison_data = {'Periodogram': {}, 'Welch': {}, 'Multitaper': {}, 'Wavelet': {}}
    psd_comparison_data['Periodogram']['psd_filt'] = psd_periodogram['filt']
    psd_comparison_data['Periodogram']['psd_u_quant'] = psd_periodogram['u_quant']
    psd_comparison_data['Periodogram']['psd_m_quant'] = psd_periodogram['m_quant']
    psd_comparison_data['Periodogram']['freqs_filt'] = psd_comparison_data['Periodogram']['freqs_u_quant'] = psd_comparison_data['Periodogram']['freqs_m_quant'] = np.fft.rfftfreq(config['signal_length'], d=1.0/config['sampling_rate'])
    win_size = 1024
    psd_comparison_data['Welch']['psd_filt'] = psd_welch['filt'][win_size]
    psd_comparison_data['Welch']['psd_u_quant'] = psd_welch['u_quant'][win_size]
    psd_comparison_data['Welch']['psd_m_quant'] = psd_welch['m_quant'][win_size]
    psd_comparison_data['Welch']['freqs_filt'] = psd_comparison_data['Welch']['freqs_u_quant'] = psd_comparison_data['Welch']['freqs_m_quant'] = np.fft.rfftfreq(win_size, d=1.0/config['sampling_rate'])
    psd_comparison_data['Multitaper']['psd_filt'] = psd_multitaper['filt']
    psd_comparison_data['Multitaper']['psd_u_quant'] = psd_multitaper['u_quant']
    psd_comparison_data['Multitaper']['psd_m_quant'] = psd_multitaper['m_quant']
    psd_comparison_data['Multitaper']['freqs_filt'] = psd_comparison_data['Multitaper']['freqs_u_quant'] = psd_comparison_data['Multitaper']['freqs_m_quant'] = np.fft.rfftfreq(config['signal_length'], d=1.0/config['sampling_rate'])

    # --- Perform Wavelet PSD calculation in Python ---
    print("Calculating Wavelet PSD in Python...")
    try:
        freqs_w_filt, psd_w_filt = estimate_psd_wavelet_py(signal_filt, config['sampling_rate'])
        freqs_w_u, psd_w_u = estimate_psd_wavelet_py(signal_u_quant, config['sampling_rate'])
        freqs_w_m, psd_w_m = estimate_psd_wavelet_py(signal_m_quant, config['sampling_rate'])
        psd_comparison_data['Wavelet']['psd_filt'] = 10 * np.log10(psd_w_filt + EPSILON)
        psd_comparison_data['Wavelet']['psd_u_quant'] = 10 * np.log10(psd_w_u + EPSILON)
        psd_comparison_data['Wavelet']['psd_m_quant'] = 10 * np.log10(psd_w_m + EPSILON)
        psd_comparison_data['Wavelet']['freqs_filt'] = freqs_w_filt
        psd_comparison_data['Wavelet']['freqs_u_quant'] = freqs_w_u
        psd_comparison_data['Wavelet']['freqs_m_quant'] = freqs_w_m
    except Exception as e:
        print(f"Could not perform Wavelet analysis: {e}")

    # --- Generate All Plots ---
    print("\n--- Generating all plots ---")
    
    # New plots as requested
    plot_time_domain_comparison(signal_orig, signal_filt, signal_u_quant, signal_m_quant, config['sampling_rate'], plot_dir)
    plot_new_acf_comparison(acf_data, config['sampling_rate'], plot_dir)
    plot_new_psd_comparison(psd_comparison_data, config['sampling_rate'], plot_dir)

    # Original plots (restored)
    plot_original_psd_comparison(psd_periodogram, config['signal_length'], config['sampling_rate'], "Periodogram", perf_data, plot_dir)
    plot_original_psd_comparison(psd_multitaper, config['signal_length'], config['sampling_rate'], "Multitaper", perf_data, plot_dir)
    plot_welch_comparison_subplots(psd_welch, perf_data, config['sampling_rate'], plot_dir)
    
    # Prepare data for distribution plot
    dist_psd_data = {
        'Periodogram': {
            'filt': psd_periodogram['filt'], 'u_quant': psd_periodogram['u_quant'], 'm_quant': psd_periodogram['m_quant']
        },
        'Welch (1024-point)': {
            'filt': psd_welch['filt'][1024], 'u_quant': psd_welch['u_quant'][1024], 'm_quant': psd_welch['m_quant'][1024]
        },
        'Multitaper': {
            'filt': psd_multitaper['filt'], 'u_quant': psd_multitaper['u_quant'], 'm_quant': psd_multitaper['m_quant']
        }
    }
    plot_full_distribution_comparison(dist_psd_data, plot_dir)
    
    # Special analysis plots
    plot_pcm_pulse_train(signal_u_quant, pcm_stream, config['sampling_rate'], plot_dir)
    plot_quantization_error(signal_filt, signal_u_quant, signal_q4, signal_q2, config['sampling_rate'], plot_dir)

    print("\nAll plots have been saved to the 'results/plots' directory.")
    print("Displaying all plots simultaneously. Close all plot windows to exit the script.")
    
    plt.show()

if __name__ == "__main__":
    main()