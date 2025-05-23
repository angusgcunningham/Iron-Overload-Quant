%Task 7: make a new file:
currentDir = pwd;

% Navigate to the data directory (assuming data is a subfolder of current directory)
dataDir = fullfile(currentDir, '1_2');

% List all files for patient 2_1 with the same sequence number
fileList = dir(fullfile(dataDir, '2_1.MR.0009.*.IMA'));

% Sort the files (they might not be in order)
fileNames = {fileList.name};
[~, sortIndex] = sort(fileNames);
fileList = fileList(sortIndex);

% Choose number of images to display
numToDisplay = min(8, length(fileList));

% Initialize arrays for images and TE values
images = cell(1, numToDisplay);
TE_values = zeros(1, numToDisplay);

% Load the selected images
for i = 1:numToDisplay
    fullPath = fullfile(fileList(i).folder, fileList(i).name);
    
    % Get DICOM info for TE value
    info = dicominfo(fullPath);
    TE_values(i) = info.EchoTime;
    
    % Read the image and convert to double
    images{i} = double(dicomread(fullPath));
    
    % Print file name and TE value
    fprintf('File %d: %s (TE = %.2f ms)\n', i, fileList(i).name, TE_values(i));
end

% Display the images in a grid
figure('Position', [100, 100, 1200, 800]);

for i = 1:numToDisplay
    subplot(2, 4, i);
    imagesc(images{i});
    colormap(gray);
    title(sprintf('TE = %.2f ms', TE_values(i)));
    axis off;
end

sgtitle('T2*-weighted MRI Images for Patient 2_1');

% Create a plot of TE vs signal intensity for a region of interest
figure;

% Define a region of interest (ROI) in the liver
% Adjust these coordinates based on where the liver appears in your images
roi_center_row = 120;  % Y-coordinate
roi_center_col = 80;  % X-coordinate
roi_radius = 8;       % Radius in pixels

% Create a circular mask for the ROI
[Y, X] = meshgrid(1:size(images{1}, 1), 1:size(images{1}, 2));
X = X';
Y = Y';
mask = ((X - roi_center_col).^2 + (Y - roi_center_row).^2) <= roi_radius^2;

% Show the ROI on the first image
figure;
imagesc(images{1});
colormap(gray);
hold on;

% Highlight the ROI
boundary = bwboundaries(mask);
plot(boundary{1}(:,2), boundary{1}(:,1), 'r', 'LineWidth', 2);
title('Selected ROI in Liver');
hold off;

% Calculate mean signal intensity in the ROI for each TE
mean_intensities = zeros(1, numToDisplay);
for i = 1:numToDisplay
    roi_pixels = images{i}(mask);
    mean_intensities(i) = mean(roi_pixels);
end

% Plot TE vs signal intensity
figure;
plot(TE_values, mean_intensities, 'bo-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Echo Time (TE) [ms]');
ylabel('Mean Signal Intensity');
title('T2* Signal Decay in Liver ROI');
grid on;

% Add an exponential fit to visualize decay
if length(TE_values) > 2
    ft = fittype('a*exp(-b*x)', 'independent', 'x');
    opts = fitoptions(ft);
    opts.StartPoint = [mean_intensities(1), 0.1];
    
    try
        fitResult = fit(TE_values', mean_intensities', ft, opts);
        
        % Calculate T2* from the fit (T2* = 1/b)
        T2star = 1/fitResult.b;
        
        % Plot the fit 
        hold on;
        x_fit = linspace(min(TE_values), max(TE_values), 100);
        y_fit = fitResult.a * exp(-fitResult.b * x_fit);
        plot(x_fit, y_fit, 'r-', 'LineWidth', 2);
        
        legend('Measured Signal', 'Exponential Fit');
        text(0.6*max(TE_values), 0.9*max(mean_intensities), ...
            sprintf('T2* = %.2f ms', T2star), 'FontSize', 12);
        hold off;
    catch
        disp('Could not fit exponential curve - need more data points');
    end
end
