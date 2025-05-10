load('pianoKeyFreqMagnitude.mat')

% File name of the .wav file
fileName = 'c3c4g4c5.wav';


% Read the .wav file
[audioData, Fs] = audioread(fileName);

% Plot settings
figure;

% Apply FFT and reduce Nyquist effect (aliasing)
n = length(audioData);         % Number of samples
y_fft = fft(audioData);        % Compute the FFT
f = (0:n-1)*(Fs/n);            % Frequency vector

% Nyquist effect reduction: limit the FFT to the positive frequencies
y_fft = y_fft(1:floor(n/2));   % Take only the first half (positive frequencies)
f = f(1:floor(n/2));           % Corresponding frequency range

% Define the maximum frequency range for the plot (0-3000 Hz)
maxFreq = 3000;

% Filter to only include frequencies within 0-3000 Hz
validFreqIdx = f <= maxFreq;   % Find frequencies <= 3000 Hz
f = f(validFreqIdx);           % Filter frequencies
y_fft = y_fft(validFreqIdx);   % Filter FFT results

% Compute magnitude and normalize by sampling rate
mag = abs(y_fft) / Fs;         % Magnitudes of the FFT divided by sampling rate
globalMax = max(mag);          % Find the global maximum magnitude

% Split frequencies into 20 Hz intervals
freqBins = 0:20:maxFreq;       % Create frequency bins of 20 Hz intervals
temp_f = [];                   % Initialize new frequency array
temp_mag = [];                 % Initialize new magnitude array

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

% Create frequency and magnitude matrix with frequencies in the first row, magnitudes in the second row
freq_magnitude = [temp_f; temp_mag];  % Create a 2-row matrix: 1st row frequencies, 2nd row magnitudes
freq_magnitude_initial = freq_magnitude;
% Plot the FFT magnitude spectrum as spikes (stem plot without circles)
stem(temp_f, temp_mag, 'Marker', 'none', 'LineWidth', 1.2); % Remove markers on spikes
title(['FFT of ' fileName ' (0-3000 Hz, Magnitude > 1/10 Global Max)']);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([0 3000]);  % Set frequency range from 0 to 3000 Hz

% Output frequency and magnitude matrix
disp('Frequency and Magnitude matrix for test.wav:');
disp(freq_magnitude);

while ~isempty(temp_f)

% Find the first note by the lowest frequency 
lowestFreq = min(temp_f);
[~, minIndex] = min(temp_f);
A4 = 440;  % Frequency of A4
semitones = round(12 * log2(lowestFreq / A4));  % Semitone difference from A4
noteNames = {'C', 'C#', 'D', 'D#', 'E', 'F', 'F#', 'G', 'G#', 'A', 'A#', 'B'};  % Musical notes

% Find the pitch corresponding to the semitone difference
octave = 4 + floor((semitones + 9) / 12);  % Calculate the octave number
noteIdx = mod(semitones + 9, 12) + 1;      % Find the note within the octave
pitch = [noteNames{noteIdx}, num2str(octave)];  % Form the pitch string (e.g., "A4")



% Retrieve the magnitude corresponding to the lowest frequency
lowestFreqMagnitude = temp_mag(minIndex);

% Display the result
disp(['pitch: ', pitch, 'magnitude: ', num2str(lowestFreqMagnitude)]);

matrixName = ['freq_mag_' pitch]; % Construct the matrix name
freq_mag_funda = eval(matrixName); % Evaluate the matrix name

% Compare the first row (frequencies) of freq_mag_funda and freq_magnitude matrices

% Get frequencies from both matrices
frequencies_funda = freq_mag_funda(1, :);  % First row of freq_mag_C3
frequencies_mag = freq_magnitude(1, :);  % First row of freq_magnitude

% Ensure that matching frequencies are found and magnitudes are updated properly

% Initialize arrays to store results
matching_frequencies = [];
matching_magnitudes_funda = [];
matching_magnitudes_mag = [];
mag_after_subtract = [];  % Array to store the updated magnitudes after subtraction

% Loop through each frequency in freq_mag_funda
for i = 1:length(frequencies_funda)
    % Find frequencies in freq_magnitude that are within 10 Hz of the current frequency in freq_mag_funda
    freq_diff = abs(frequencies_mag - frequencies_funda(i));
    
    % If any frequency matches within 10 Hz
    if any(freq_diff <= 5)
        % Get the index of the closest matching frequency
        [~, matchIdx] = min(freq_diff);
        
        % Store the matching frequency and its corresponding magnitudes
        matching_frequencies = [matching_frequencies, frequencies_funda(i)];
        matching_magnitudes_funda = [matching_magnitudes_funda, freq_mag_funda(2, i)];  % Magnitude from freq_mag_funda
        matching_magnitudes_mag = [matching_magnitudes_mag, freq_magnitude(2, matchIdx)];  % Magnitude from freq_magnitude
        
        % Subtract the magnitudes and store the result
        correspondingMagnitude = matching_magnitudes_mag(end); % Ensure this value is set correctly
        subtracted_mag = correspondingMagnitude - lowestFreqMagnitude*matching_magnitudes_funda(end);
        mag_after_subtract = [mag_after_subtract, subtracted_mag];
        
        % Update temp_mag with the new subtracted magnitude
        temp_mag(matchIdx) = subtracted_mag; % Ensure matchIdx is valid
    end
end


% Create a matrix of matching frequencies and their corresponding magnitudes from both sources
matching_results = [matching_frequencies; matching_magnitudes_funda; matching_magnitudes_mag; mag_after_subtract];

% Output the matching frequencies and magnitudes
disp('Matching frequencies and their magnitudes from freq_mag_funda and freq_magnitude:');
disp(matching_results);



% Display the updated temp_mag array
disp('Updated temp_mag after subtraction:');
disp(temp_mag);


% Remove frequencies with magnitudes less than 0.01
threshold = 0.1; % Magnitude threshold
validIndices = temp_mag >= threshold; % Find indices where magnitude is above the threshold

% Update temp_f and temp_mag by keeping only valid indices
temp_f = temp_f(validIndices);
temp_mag = temp_mag(validIndices);

% Plot the updated FFT magnitude spectrum
figure;
stem(temp_f, temp_mag, 'Marker', 'none', 'LineWidth', 1.2); % Remove markers on spikes
title(['Updated FFT of ' fileName ' (0-3000 Hz, Magnitude > 0.01)']);
xlabel('Frequency (Hz)');
ylabel('Normalized Magnitude');
xlim([0 3000]);  % Set frequency range from 0 to 3000 Hz
ylim([0 1]);

% Output the updated frequency and magnitude matrix
freq_magnitude = [temp_f; temp_mag]; % Create a 2-row matrix: 1st row frequencies, 2nd row magnitudes
disp('Updated Frequency and Magnitude matrix:');
disp(freq_magnitude);
end

disp('Notes and their corresponding LowestFreqMagnitude:');
for i = 1:size(notesWithMagnitudes, 1)
    disp([notesWithMagnitudes{i, 1}, ' with magnitude: ', num2str(notesWithMagnitudes{i, 2})]);
end