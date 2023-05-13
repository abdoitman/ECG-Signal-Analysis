# ECG Signal Analysis
Used MATLAB to clean an ECG signal and extract useful information about it.

## Main Points
  1. Removed unwanted low frequncy noise caused by the movement of the patient during the measurement by **applying a high pass filter**.
  2. ECG is contaminated with 50 Hz noise from the power outlet, which overlaps the ECG. This interference was removed by **applying a notch filter at 50 Hz**.
  3. To furthur improve the signal to noise ratio, **various low pass filters were applied with different cutoff frequencies** to determine the best.
  4. Heart rate was calculated using **ACF (autocorrelation function) of the function**.
  5. Implemented PAn-Tompkins algorithm to detect QRS complex of the signal.
