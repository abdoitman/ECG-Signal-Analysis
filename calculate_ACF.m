function heart_rate_acf = calculate_ACF(ECGname, ecg)
    samples = length(ecg);
    Fs = 500;

    acf = xcorr(ecg);
    Shift_vector = linspace(-samples, samples, samples*2-1);
    figure;
    plot(Shift_vector, acf)
    title("Auto correlation function for " + ECGname)
    xlabel("Shift")
    ylabel("ACF")
    [max_peaks, max_locs] = findpeaks(acf);
    [global_max, global_loc] = max(acf);
    local_max = max(max_peaks(max_locs > global_loc));
    local_loc = max_locs(max_peaks == local_max & max_locs > global_loc);
    heart_rate_acf = Fs / (local_loc - global_loc)*60;
end