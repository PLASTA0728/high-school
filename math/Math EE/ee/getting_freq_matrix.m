% List of piano keys from A0 to C8 (88 keys)
% Assuming the .wav files are named accordingly, e.g., 'A0.wav', 'Bb0.wav', ..., 'C8.wav'

pianoKeys = {'A0', 'Bb0', 'B0', 'C1', 'Db1', 'D1', 'Eb1', 'E1', 'F1', 'Gb1', 'G1', 'Ab1', ...
             'A1', 'Bb1', 'B1', 'C2', 'Db2', 'D2', 'Eb2', 'E2', 'F2', 'Gb2', 'G2', 'Ab2', ...
             'A2', 'Bb2', 'B2', 'C3', 'Db3', 'D3', 'Eb3', 'E3', 'F3', 'Gb3', 'G3', 'Ab3', ...
             'A3', 'Bb3', 'B3', 'C4', 'Db4', 'D4', 'Eb4', 'E4', 'F4', 'Gb4', 'G4', 'Ab4', ...
             'A4', 'Bb4', 'B4', 'C5', 'Db5', 'D5', 'Eb5', 'E5', 'F5', 'Gb5', 'G5', 'Ab5', ...
             'A5', 'Bb5', 'B5', 'C6', 'Db6', 'D6', 'Eb6', 'E6', 'F6', 'Gb6', 'G6', 'Ab6', ...
             'A6', 'Bb6', 'B6', 'C7', 'Db7', 'D7', 'Eb7', 'E7', 'F7', 'Gb7', 'G7', 'Ab7', ...
             'A7', 'Bb7', 'B7', 'C8'};

% Initialize storage for audio data, sampling rates, frequencies, and magnitudes
audioData = cell(1, 88); % Store the audio data
Fs = zeros(1, 88);       % Store the sampling rates
freq_magnitude = cell(1, 88); % Store frequency and magnitude matrices for each file

% Define the maximum frequency range for processing (0-3000 Hz)
maxFreq = 3000;

% Loop to read all 88 .wav files
for i = 1:88
    % Read the .wav file
    fileName = [pianoKeys{i}, '.mp3']; % Construct the file name
    [audioData{i}, Fs(i)] = audioread(fileName);
end

freq_mag_map = containers.Map;

% Loop to process each file and compute FFT
for i = 1:88
    % Apply FFT and reduce Nyquist effect (aliasing)
    n = length(audioData{i});         % Number of samples
    y_fft = fft(audioData{i});        % Compute the FFT
    f = (0:n-1)*(Fs(i)/n);            % Frequency vector

    % Nyquist effect reduction: limit the FFT to the positive frequencies
    y_fft = y_fft(1:floor(n/2));      % Take only the first half (positive frequencies)
    f = f(1:floor(n/2));              % Corresponding frequency range

    % Filter to only include frequencies within 0-3000 Hz
    validFreqIdx = f <= maxFreq;      % Find frequencies <= 3000 Hz
    f = f(validFreqIdx);              % Filter frequencies
    y_fft = y_fft(validFreqIdx);      % Filter FFT results
    
    % Compute magnitude and normalize by sampling rate
    mag = abs(y_fft) / Fs(i);         % Magnitudes of the FFT divided by sampling rate
    globalMax = max(mag);             % Find the global maximum magnitude

    % Split frequencies into 20 Hz intervals
    freqBins = 0:20:maxFreq;          % Create frequency bins of 20 Hz intervals
    temp_f = [];                      % Initialize new frequency array
    temp_mag = [];                    % Initialize new magnitude array
    
    % Loop over each 20 Hz interval
    for j = 1:length(freqBins)-1
        % Find the indices of frequencies within the current 20 Hz bin
        binIdx = (f >= freqBins(j)) & (f < freqBins(j+1));
        bin_frequencies = f(binIdx);  % Frequencies in the current bin
        bin_magnitudes = mag(binIdx); % Magnitudes in the current bin
        
        % If there are any frequencies in this bin, find the maximum magnitude
        if ~isempty(bin_magnitudes)
            [maxMag, maxIdx] = max(bin_magnitudes); % Find max magnitude in the bin
            
            % Check if the max magnitude is greater than 1/10th of the global maximum
            if maxMag >= globalMax / 10
                % Add the maximum frequency and its magnitude to the new arrays
                temp_f = [temp_f, bin_frequencies(maxIdx)];  % Maximum frequency in this bin
                temp_mag = [temp_mag, maxMag / globalMax];   % Normalized maximum magnitude in this bin
            end
        end
    end
    
    % Check minimum frequency and cancel spikes based on neighboring ranges
    minFreq = min(temp_f); % Find the minimum frequency in the current results
    if minFreq > 40
        % Identify and remove lower magnitude spikes if frequency difference < 40 Hz
        for j = 1:length(temp_f)-1
            if abs(temp_f(j+1) - temp_f(j)) < 40
                if temp_mag(j) < temp_mag(j+1)
                    temp_f(j) = NaN;  % Remove lower frequency
                    temp_mag(j) = NaN; % Remove lower magnitude
                else
                    temp_f(j+1) = NaN; % Remove lower frequency
                    temp_mag(j+1) = NaN; % Remove lower magnitude
                end
            end
        end
        
        % Remove NaN values (lower spikes)
        temp_f = temp_f(~isnan(temp_f));
        temp_mag = temp_mag(~isnan(temp_mag));
    end
    
    % Normalize magnitudes by the magnitude of the first spike
    if ~isempty(temp_mag)
        firstSpikeMagnitude = temp_mag(1); % Use the first spike magnitude for normalization
        temp_mag = temp_mag / firstSpikeMagnitude; % Normalize all magnitudes
    end
    
dynamic_field = sprintf('freq_mag_%s', pianoKeys{i}); % Generate field name
    freq_mag_struct.(dynamic_field) = freq_magnitude{i};  % Store the magnitude
end


