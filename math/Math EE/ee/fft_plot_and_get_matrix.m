% File names of the four specific .wav files
fileNames = {'C3.wav','C4.wav', 'C5.wav', 'G4.wav'};

% Initialize storage for audio data, sampling rates, frequencies, and magnitudes
audioData = cell(1,4); % Store the audio data
Fs = zeros(1,4);       % Store the sampling rates
freq_magnitude = cell(1,4); % Store frequency and magnitude matrices for each file

% Define the maximum frequency range for the plots (0-3000 Hz)
maxFreq = 3000;

% Loop to read the four .wav files
for i = 1:4
    % Read the .wav file
    [audioData{i}, Fs(i)] = audioread(fileNames{i});
end

% Plot settings
figure;

% Loop to process each file and compute FFT
for i = 1:4
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
    new_f = [];                       % Initialize new frequency array
    new_mag = [];                     % Initialize new magnitude array
    
    % Storage for intermediate results
    temp_f = [];  % Store frequencies
    temp_mag = []; % Store magnitudes
    
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
        % Create a copy of the temporary arrays for modification
        final_f = temp_f;
        final_mag = temp_mag;
        
        % Identify and remove lower magnitude spikes if frequency difference < 40 Hz
        for j = 1:length(final_f)-1
            if abs(final_f(j+1) - final_f(j)) < 40
                % Compare magnitudes and remove the lower one
                if final_mag(j) < final_mag(j+1)
                    final_f(j) = NaN;  % Mark the lower frequency to be removed
                    final_mag(j) = NaN; % Mark the lower magnitude to be removed
                else
                    final_f(j+1) = NaN; % Mark the lower frequency to be removed
                    final_mag(j+1) = NaN; % Mark the lower magnitude to be removed
                end
            end
        end
        
        % Remove NaN values (lower spikes)
        final_f = final_f(~isnan(final_f));
        final_mag = final_mag(~isnan(final_mag));
        
        % Update the results
        temp_f = final_f;
        temp_mag = final_mag;
    end
    
    % Normalize magnitudes by the magnitude of the first spike
    if ~isempty(temp_mag)
        firstSpikeMagnitude = temp_mag(1); % Use the first spike magnitude for normalization
        temp_mag = temp_mag / firstSpikeMagnitude; % Normalize all magnitudes
    end
    
    % Create frequency and magnitude matrix with frequencies in first row, magnitudes in second row
    freq_magnitude{i} = [temp_f; temp_mag];  % Create a 2-row matrix: 1st row frequencies, 2nd row magnitudes

    % Plot the FFT magnitude spectrum as spikes (stem plot without circles)
    subplot(2, 2, i); % Arrange the four plots in a 2x2 grid
    stem(temp_f, temp_mag, 'Marker', 'none', 'LineWidth', 1.2); % Remove markers on spikes
    title(['FFT of ' fileNames{i} ' (0-3000 Hz, Normalized Magnitude)']);
    xlabel('Frequency (Hz)');
    ylabel('Normalized Magnitude');
    xlim([0 3000]);  % Set frequency range from 0 to 3000 Hz
end

% Output frequency and magnitude matrices
freq_mag_C3 = freq_magnitude{1};  % Frequencies and magnitudes for C3.wav
freq_mag_C4 = freq_magnitude{2};  % Frequencies and magnitudes for C4.wav
freq_mag_C5 = freq_magnitude{3};  % Frequencies and magnitudes for C5.wav
freq_mag_G4 = freq_magnitude{4}; % Frequencies and magnitudes for chord.wav

% Display the frequency and magnitude matrices in the command window (optional)
disp('C3.wav:');
disp(freq_mag_C3);
disp('C4.wav:');
disp(freq_mag_C4);
disp('C5.wav:');
disp(freq_mag_C5);
disp('G4.wav:');
disp(freq_mag_G4);
