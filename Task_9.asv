%% Task 9 Process All Patients with 4 Fixed ROIs (each ~200 pixels)

% Define current directory and patient folder names
currentDir = pwd;
patientFolders = {'1_2', '2_7', '3_8', '4_15'};
numPatients = length(patientFolders);

% Preallocate cells to store TE, and estimated parameter maps for each patient
TE_all = cell(numPatients, 1);
a_all = cell(numPatients, 1);  % Estimated a (initial signal)
r_all = cell(numPatients, 1);  % Estimated r (T2* computed as 1./r)

% Use a file pattern for DICOM images
filePattern = '*.IMA';

for p = 1:numPatients
    % Define the data directory for the patient
    dataDir = fullfile(currentDir, patientFolders{p});
    
    % List and sort all DICOM files in the folder
    fileList = dir(fullfile(dataDir, filePattern));
    fileNames = {fileList.name};
    [~, sortIndex] = sort(fileNames);
    fileList = fileList(sortIndex);
    
    % Number of echoes/images
    numEchoes = length(fileList);
    
    % Preallocate arrays to store TE values and images
    TE_values = zeros(1, numEchoes);
    images_patient = cell(1, numEchoes);
    
    % Load each DICOM file, extract TE, and read the image
    for i = 1:numEchoes
        fullPath = fullfile(fileList(i).folder, fileList(i).name);
        info = dicominfo(fullPath);
        TE_values(i) = info.EchoTime;  % TE in ms
        images_patient{i} = double(dicomread(fullPath));
    end
    TE_all{p} = TE_values;
    
    % Get image dimensions (assume all images have same dimensions)
    [Nrow, Ncol] = size(images_patient{1});
    bands = numEchoes;
    
    % Stack images into a 3D volume: size [Nrow x Ncol x bands]
    vol = zeros(Nrow, Ncol, bands);
    for i = 1:bands
        vol(:,:,i) = images_patient{i};
    end
    
    % Reshape the volume to create a 2D data matrix for estimation:
    % Each column is the time series for one pixel.
    yReshaped = reshape(vol, Nrow*Ncol, bands)';  % dimensions: [bands x (Nrow*Ncol)]
    
    % Set TV regularization parameters (adjust if needed)
    lambdaA = 1e-5;
    lambdaR = 1e-5;
    
    % Run the ADMM-based T2* estimation
    [a_est, r_est] = relaxationEst(yReshaped, TE_values, Nrow, Ncol, lambdaA, lambdaR);
    
    % Reshape estimated vectors back into 2D parameter maps and store them.
    a_all{p} = reshape(a_est, Nrow, Ncol);
    r_all{p} = reshape(r_est, Nrow, Ncol);  % T2* will be computed as 1./r

    % Enforce a lower bound so r never goes zero or negative:
    eps_floor = 1e-6;  
    r_all{p}(r_all{p} < eps_floor) = eps_floor;
    % ————————————————

end

%% ROI Analysis: Define 4 Fixed ROIs (each ~200 pixels) around a reference ROI

% Define the reference ROI center (from previous Task 7) and an offset to place 4 ROIs
ref_row = 120;   % Reference row (adjust based on your image)
ref_col = 80;   % Reference column
offset = 8;     % Distance (in pixels) to shift ROI centers relative to the reference

% Define 4 ROI centers (for example, in the four quadrants around the reference center)
roi_centers = [ref_row - offset, ref_col - offset;  % top-left
               ref_row - offset, ref_col + offset;  % top-right
               ref_row + offset, ref_col - offset;  % bottom-left
               ref_row + offset, ref_col + offset]; % bottom-right
           
% The area of each ROI is 200 pixels so the radius of a circular ROI is:
roi_radius = sqrt(200/pi);

% Preallocate matrices to store mean values for each ROI and each patient
mean_a_ROIs = zeros(numPatients, 4);   % Estimated a in 4 ROIs
mean_T2_ROIs = zeros(numPatients, 4);    % Estimated T2* (1./r)

for p = 1:numPatients
    % Get the current patient's estimated parameter maps
    a_img = a_all{p};
    T2_img = 1 ./ r_all{p};  % Compute T2* map
    
    % Get image dimensions and create a grid for mask creation
    [rows, cols] = size(a_img);
    [X, Y] = meshgrid(1:cols, 1:rows);
    
    % For each of the 4 ROIs, create a circular mask and compute the mean values
    for roi = 1:4
         center_row = roi_centers(roi, 1);
         center_col = roi_centers(roi, 2);
         mask = ((X - center_col).^2 + (Y - center_row).^2) <= roi_radius^2;
         mean_a_ROIs(p, roi) = mean(a_img(mask));
         mean_T2_ROIs(p, roi) = mean(T2_img(mask));
    end
    
    % Optionally, display the T2* map with ROI boundaries overlaid
    figure;
    imagesc(T2_img);
    colormap(hsv);
    colorbar;
    axis image off;
    hold on;
    for roi = 1:4
         center_row = roi_centers(roi, 1);
         center_col = roi_centers(roi, 2);
         mask = ((X - center_col).^2 + (Y - center_row).^2) <= roi_radius^2;
         boundaries = bwboundaries(mask);
         for k = 1:length(boundaries)
             plot(boundaries{k}(:,2), boundaries{k}(:,1), 'r', 'LineWidth', 2);
         end
    end
    title(sprintf('Estimated T2* with 4 ROI Boundaries (Patient %s)', patientFolders{p}));
    hold off;
end

%% Display Parameter Maps for All Patients (with clamped T2* display)
figure('Position', [100, 100, 1200, 800]);
for p = 1:numPatients
    % --- Amplitude map ---
    subplot(numPatients, 2, (p-1)*2+1);
    imagesc(a_all{p});
    colormap(gray);
    colorbar;
    axis image off;
    title(sprintf('Estimated a for Patient %s', patientFolders{p}));
    
    % --- T2* map ---
    subplot(numPatients, 2, (p-1)*2+2);
    % Clamp extremes
    T2_map = 1 ./ r_all{p};
    T2_map(T2_map < 0)  = NaN;
    T2_map(T2_map > 30) = 30;
    imagesc(T2_map, [0 30]);
    colormap(hot);
    h = colorbar;
    h.Ticks = 0:5:30;
    axis image off;
    
    % Compute the mean T2* across your four fixed ROIs
    meanT2 = mean(mean_T2_ROIs(p,:));
    
    % Two-line title matching Task 11 style
    title({
      sprintf('Patient %s - T2* Map (ADMM)', patientFolders{p}), ...
      sprintf('Mean: %.2f ms', meanT2)
    });
end



%% Display ROI Mean Values as Tables

roi_labels = {'ROI1','ROI2','ROI3','ROI4'};
T2_Table = array2table(mean_T2_ROIs, 'VariableNames', roi_labels, 'RowNames', patientFolders);
a_Table = array2table(mean_a_ROIs, 'VariableNames', roi_labels, 'RowNames', patientFolders);

disp('Mean Estimated T2* Values in 4 Fixed ROIs:');
disp(T2_Table);

disp('Mean Estimated a Values in 4 Fixed ROIs:');
disp(a_Table);

%% --- Display the image and overlay the 4 ROIs ---
figure;
imshow(I,[]);
colormap(gray);
hold on;
for k=1:4
    ctr = roi_centers(k,:);
    viscircles([ctr(2), ctr(1)], roi_radius, 'Color','r','LineWidth',1);
end
title('Patient 2\_7: Four ROIs overlaid on first echo');
hold off;