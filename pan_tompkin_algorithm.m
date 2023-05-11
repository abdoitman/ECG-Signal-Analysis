function [QRS_duration,heart_rate]=pan_tompkin_algorithm(ECGname, ecg,fs, T_vector)
    %%%% Constants %%%%
    ecg_len = length(ecg);
    Fs = fs;

    %%%% Low pass filter %%%%
    b_lpf = [1 0 0 0 0 0 -2 0 0 0 0 0 1];
    a_lpf = [1 -2 1] * 32;
    
    lpf = dfilt.df2(b_lpf, a_lpf);
    
    ecg_filtered_2 = filter(lpf, ecg);
    ecg_filtered_lpf = ecg_filtered_2 - mean(ecg_filtered_2);
%     ecg_filtered_lpf = ecg_filtered_2 / max(abs(ecg_filtered_2));

    %%%% High pass filter %%%%
    b_hpf = [-1, zeros(1,15), 32, -32, zeros(1,14), 1];
    a_hpf = [1, -1] * 32;

    bpf = dfilt.df2(b_hpf, a_hpf);

    ecg_filtered_2 = filter(bpf, ecg_filtered_lpf);
    ecg_filtered_bpf = [ecg_filtered_2(1:40) *0.25;
        ecg_filtered_2(41:end)];
%     ecg_filtered_bpf = ecg_filtered_2 / max(abs(ecg_filtered_2));

    figure;
    subplot(3,2,[1 2])
    plot(T_vector, ecg)
    title("Original ECG for " + ECGname + "signal")
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")

    subplot(3,2,3)
    plot(T_vector, ecg_filtered_bpf)
    title("ECG signal after applying BPF")
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")

    %%%% Diffrintiator ####
    b_d = [2 1 0 -1 2];
    a_d = [1] * 8;

    hd = dfilt.df2(b_d, a_d);

    ecg_diff = filter(hd, ecg_filtered_bpf);
%     ecg_diff = ecg_diff / max(abs(ecg_diff));

    subplot(3,2,4)
    plot(T_vector, ecg_diff)
    title("ECG signal after the differentiation operation")
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")

    %%%% Squaring %%%%
    ecg_squared = ecg_diff .^ 2;
%     ecg_squared = ecg_squared / max(abs(ecg_squared));

    subplot(3,2,5)
    plot(T_vector, ecg_squared)
    title("ECG signal after the squaring operation")
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")

    %%%% Integrating Window %%%%
    ecg_integ_pad = [zeros(1,29) ecg_squared' zeros(1,29)];
    for i = 30:length(ecg_integ_pad) - 29
        ecg_integ(i - 29) = sum(ecg_integ_pad(i-29:i))/30;
    end
    ecg_integ = ecg_integ';
%     ecg_integ = ecg_integ/max(abs(ecg_integ));

    subplot(3,2,6)
    plot(T_vector, ecg_integ)
    title("ECG signal after the squaring operation")
    xlabel("Time (seconds)")
    ylabel("Voltage (V)")

    %%% Thresholding & Detecting RQS %%%
    TH = mean(ecg_integ);
    ecg_th = zeros(ecg_len, 1);
    w = (ecg_integ>(TH));
    ecg_th(w) = 1;
    x = find(diff([0 w'])== 1);
    y = find(diff([w' 0])== -1);
    x = x - 55;
    y = y - 55;

    for i= 1:length(x)
        [R_val(i), R_loc(i)] = max(ecg(x(i):y(i)));
        R_loc(i) = R_loc(i) - 1 + x(i);
        [Q_val(i), Q_loc(i)] = min(ecg(R_loc(i):-1:R_loc(i)-18 ));
        Q_loc(i) = R_loc(i) - Q_loc(i) +1;
        [S_val(i), S_loc(i)] = min(ecg(R_loc(i):R_loc(i)+15));
        S_loc(i) = R_loc(i) + S_loc(i) -1;
    end

   %%%% R Peak localization %%%%
   figure;
   subplot(2,1,1);
   plot(T_vector, ecg, T_vector(R_loc), R_val, 'r^', T_vector(S_loc), S_val, '*', T_vector(Q_loc), Q_val, 'o');
   title("RQS detection in one period of " + ECGname + " signal")
   xlim([0.5 2.5]);
   legend("ECG", "R", "S", "Q");
   
   subplot(2,1,2);
   plot(T_vector, ecg, T_vector(R_loc), R_val, 'r^', T_vector(S_loc), S_val, '*', T_vector(Q_loc), Q_val, 'o');
   title("RQS detection in " + ECGname + " signal")
   legend("ECG", "R", "S", "Q");

   %%% Calculation of QRS duration and heart rate %%%
   count = -1;
   for i = 1:ecg_len-1
       if(ecg_th(i) == 0 && ecg_th(i+1) >= 1 && ecg_th(i+2) >= 1)
           count = count + 1;
       end
   end

   ecg_pbm = diff([ecg_th; 0]);
   x = find(ecg_pbm > 0);
   y = find(ecg_pbm < 0);
   z = y - x;

   QRS_duration = mean(z) * (1/Fs);
   heart_rate = count * (1/Fs) * (ecg_len); 

   disp("QRS Duration = " + QRS_duration)
    disp("Heart Rate from Pan Tompkins = " + heart_rate)
end