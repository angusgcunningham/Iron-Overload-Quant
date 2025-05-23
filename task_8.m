%% Task 8: Plot TE vs Signal Intensity for Patients 1_2 and 3_8

% Define current directory (adjust if needed)
currentDir = pwd;

% Define folders for two patients – update these folder names if necessary
patientFolders = {'1_2', '3_8'};
numPatients = length(patientFolders);

% Preallocate cell arrays to store TE values and mean intensities for each patient
TE_all = cell(numPatients,1);
meanIntensities_all = cell(numPatients,1);

% Define the file pattern; adjust this if needed for your DICOM files
filePattern = '*.IMA';

for p = 1:numPatients
    % Set the data directory for the current patient
    dataDir = fullfile(currentDir, patientFolders{p});
    
    % List all DICOM files in the patient folder
    fileList = dir(fullfile(dataDir, filePattern));
    
    % Sort the file names (to keep TE order consistent)
    fileNames = {fileList.name};
    [~, sortIndex] = sort(fileNames);
    fileList = fileList(sortIndex);
    
    % Choose number of images (echoes) to use (for example, 8 echoes)
    numToDisplay = min(8, length(fileList));
    
    % Initialize arrays to store images and TE values for current patient
    TE_values = zeros(1, numToDisplay);
    images_patient = cell(1, numToDisplay);
    
    % Load images and TE values
    for i = 1:numToDisplay
        fullPath = fullfile(fileList(i).folder, fileList(i).name);
        info = dicominfo(fullPath);
        TE_values(i) = info.EchoTime; % TE in ms
        images_patient{i} = double(dicomread(fullPath));
    end
    
    % Create the ROI mask based on the size of the images (assuming all are same size)
    [rows, cols] = size(images_patient{1});
    [X, Y] = meshgrid(1:cols, 1:rows); % X: columns, Y: rows
    mask = ((X - roi_center_col).^2 + (Y - roi_center_row).^2) <= roi_radius^2;
    
    % Optionally, visualize the ROI on the first image for confirmation:
    figure;
    imagesc(images_patient{1});
    colormap(gray);
    hold on;
    boundary = bwboundaries(mask);
    plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', 2);
    title(sprintf('ROI for patient %s', patientFolders{p}));
    hold off;
    
    % Compute the mean intensity within the ROI for each echo (TE)
    mean_intensities = zeros(1, numToDisplay);
    for i = 1:numToDisplay
        currentImage = images_patient{i};
        mean_intensities(i) = mean(currentImage(mask));
    end
    
    % Save TE values and computed mean intensities for the current patient
    TE_all{p} = TE_values;
    meanIntensities_all{p} = mean_intensities;
    
    % Optionally, you can print some diagnostic info
    fprintf('Patient %s:\n', patientFolders{p});
    disp(table(TE_values', mean_intensities', 'VariableNames', {'TE_ms','MeanIntensity'}));
end

%% Plot TE vs. Mean Signal Intensity for Both Patients on the Same Graph

figure;
markers = {'bo-', 'rs-'};
hold on;
for p = 1:numPatients
    plot(TE_all{p}, meanIntensities_all{p}, markers{p}, 'LineWidth',2, 'MarkerSize',8);
end
hold off;
xlabel('Echo Time (TE) [ms]');
ylabel('Mean Signal Intensity');
title('TE vs. Mean Signal Intensity for Patients 1_2 and 3_8');
legend(patientFolders, 'Location','best');
grid on;
