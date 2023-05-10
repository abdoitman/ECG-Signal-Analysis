%% Constants
load("ecg.mat")
Fs  = 500;
T = 1 / Fs;
samples = length(ecg);
samples_per_Hz = samples / Fs;
EKG1 = 500;
ecg = ecg / EKG1;

%% Q1
T_vector = linspace(0, samples * T, samples);
figure;
plot(T_vector, ecg);
title("ECG Signal")
xlabel("Time (seconds)");
ylabel("Voltage (V)");


ecg1_unfiltered = fftshift(fft(ecg));
F_vector = linspace(-Fs/2, Fs/2, 4170);

H = ones(4170,1);
H([samples/2 - floor(0.5 * samples_per_Hz):samples/2 + ceil(0.5*samples_per_Hz)]) = 0;
% d = linspace(ceil(samples/2 - 0.5 * samples_per_Hz), floor(samples/2 + 0.5 * samples_per_Hz));
% H(d) = 0;

ecg1_filtered = ecg1_unfiltered .* H;
figure;
plot(F_vector, abs(ecg1_filtered))
title("Filtered ECG in the frequency domain")
xlabel("Frequency (Hz)")
ylabel("ECG signal in frequency domain")
figure;
plot(F_vector, H);
title("Transfer function in the frequency domain")
xlabel("Frequency (Hz)")
ylabel("Transfer function in frequency domain")

ecg1 = real(ifft(ifftshift((ecg1_filtered))));
H_in_time = abs(ifft(ifftshift(H)));
figure;
plot(T_vector, ecg1)
title("Filtered ECG in the time domain")
xlabel("Time (seconds)")
ylabel("Voltage (V)")
figure;
plot(T_vector, H_in_time);
title("Impulse response of transfer function")
xlabel("Time")
ylabel("Impulse response")

%% Q2

f0 = 50; %Hz
wo = f0/(Fs/2);
bw = wo/35;

[b, a] = iirnotch(wo, bw);
ecg2 = filter(b, a , ecg1);

figure;
plot(T_vector, ecg2)
title('ECG signal after removing the noise at 50 Hz')
xlabel('Time (seconds)')
ylabel('Voltage (V)')

figure;
plot(F_vector, abs(fftshift(fft(ecg2))))

%% Q3

Cutoff_frequencies = [60, 50, 40, 30, 20, 10];
T_one_period = [ceil(0.5*Fs):ceil(Fs*1.5)];

ecg3_freq = fftshift(fft(ecg2));

figure;
for i = 1:length(Cutoff_frequencies)
    ecg3_freq_filtered = ecg3_freq;
    frequency_point = Cutoff_frequencies(i) * 4170 / 500;
    ecg3_freq_filtered([1:(samples/2) - ceil(frequency_point) (samples/2) + floor(frequency_point):end]) = 0;
    ecg3 = real(ifft(ifftshift(ecg3_freq_filtered)));

    subplot(2,3,i)
    plot(T_vector(T_one_period), ecg3(T_one_period));
    hold on
    plot(T_vector(T_one_period), ecg2(T_one_period));
    legend(["Filtered Signal", "Original Signal"])
    title("One period of the filtered signal vs the original at cutoff frequence = " + Cutoff_frequencies(i))
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")
    hold on
end
hold off

figure;
for i = 1:length(Cutoff_frequencies)
    ecg3_freq_filtered = ecg3_freq;
    frequency_point = Cutoff_frequencies(i) * 4170 / 500;
    ecg3_freq_filtered([1:(samples/2) - ceil(frequency_point) (samples/2) + floor(frequency_point):end]) = 0;

    subplot(2,3,i)
    plot(F_vector, abs(ecg3_freq_filtered));
    title("Filtered signal with cutoff frequncy = " + Cutoff_frequencies(i))
    xlabel("Frequency (Hz)")
end

% Best cutoff frequency is at 40 Hz

%% Q4

acf = xcorr(ecg3);
[max_peaks, max_locs] = findpeaks(acf);
[global_max, global_loc] = max(acf);
local_max = max(max_peaks(max_locs > global_loc));
local_loc = max_locs(max_peaks == local_max & max_locs > global_loc);
heart_rate = Fs / (local_loc - global_loc)*60;

%% Q5

