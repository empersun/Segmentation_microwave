%%%%%%%5 EXAMPLE FROM GITHUB %%%%%%%%%%%

%% Load sample data_high_loc3_no_refl (antenna locations, frequencies and signals)
clc;
frequencies = dlmread('C:\MCGILL\kunwei_code\Microwave\MERIT-master\data_high_loc3\frequencies4.csv');
antenna_locations = dlmread('C:\MCGILL\kunwei_code\Microwave\MERIT-master\data_high_loc3\antenna_locations.csv');
channel_names = dlmread('C:\MCGILL\kunwei_code\Microwave\MERIT-master\data_high_loc3\channel_names.csv');

load('data_high_loc3\scans.mat')
% scan1 = dlmread('data_high_loc3_no_refl/B0_P3_p000.csv');
% scan2 = dlmread('data_high_loc3_no_refl/B0_P3_p036.csv');
%fre1
% scan1=scan11(1:81,:);
% scan2=scan22(1:81,:);
% fre2
% scan1=scan11(1:41,:);
% scan2=scan22(1:41,:);
% fre3
% scan1=scan11(1:61,:);
% scan2=scan22(1:61,:);
% fre4
scan1=scan1(22:61,:);
scan2=scan2(22:61,:);
% 
% %% Perform rotation subtraction
signals = scan2 - scan1;
% signals2 = scan2_freq2 - scan1_freq2;
% signals3 = scan2_freq3 - scan1_freq3;
% signals4 = scan2_freq4 - scan1_freq4;

% % %% Calculate delays for synthetic focusing
% % delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
% %   'relative_permittivity', 6);
% % 
% % %% Perform imaging
% % img = abs(merit.beamform(signals, frequencies, points, delays, ...
% %         merit.beamformers.DAS));

% % % Generate imaging domain
% [points, axes_] = merit.domain.hemisphere('radius', 0.075, 'resolution', 2e-3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define parameters for the cylindrical domain
radius = 0.075;  % Radius of the cylinder
resolution = 2e-3;  % Resolution
height = 0.08;  % Height of the cylinder (adjust as needed)

% Define the constant z value for the slice
z_slice_value = 0.05;

% Create a grid of points within the cylindrical domain
[X, Y, Z] = meshgrid(-radius:resolution:radius, -radius:resolution:radius, 0:resolution:height);

% Filter points to be within the radius
mask = sqrt(X.^2 + Y.^2) <= radius;
X = X(mask);
Y = Y(mask);
Z = Z(mask);

% Combine the points into a single array
points = [X(:), Y(:), Z(:)];

%% Calculate delays for synthetic focusing
delays = merit.beamform.get_delays(channel_names, antenna_locations, ...
    'relative_permittivity', 9);

%% Perform imaging
img = abs(merit.beamform(signals, frequencies, points, delays, ...
    merit.beamformers.DAS));

% % Define the constant z value for the slice
% z_slice_value = 0.05;

% Find the indices of the points at the desired z value
slice_indices = abs(Z - z_slice_value) < resolution/2;
slice_points = points(slice_indices, :);

% Get the corresponding axes for the slice
slice_axes_x = unique(slice_points(:, 1));
slice_axes_y = unique(slice_points(:, 2));

% Extract the slice of the image
im_slice = merit.visualize.get_slice(img, slice_points, {slice_axes_x, slice_axes_y}, 'z', z_slice_value);

% Display the image
imagesc(slice_axes_x, slice_axes_y, im_slice');
xlabel('X-axis');
ylabel('Y-axis');
title(['Slice at Z = ', num2str(z_slice_value)]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% KUNWEI %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% Extract the x, y, and z coordinates
x_locations = antenna_locations(:, 1);
y_locations = antenna_locations(:, 2);
z_locations = antenna_locations(:, 3); % z_locations are not used for 2D plotting

% % % Get the slice of the image
% % im_slice = merit.visualize.get_slice(img, points, axes_, 'z', 0.05);
% % 
% % % Display the image
% % imagesc(axes_{1:2}, im_slice');

% % Set the axes to have the same length for data units along each axis and fit tightly around the data
% axis image;
% 
% % Get the current axes handle
% ax = gca;
% 
% % Determine the range of values for the ticks based on the current axis limits
% x_limits = get(ax, 'XLim');
% y_limits = get(ax, 'YLim');
% min_limit = min(x_limits(1), y_limits(1));
% max_limit = max(x_limits(2), y_limits(2));
% 
% % Create a set of tick marks
% ticks = linspace(min_limit, max_limit, 5);
% 
% % Set the same tick marks for both the x and y axes
% set(ax, 'XTick', ticks, 'YTick', ticks);
% 
% % Hold the current plot to overlay the antenna locations
% hold on;
% 
% % Plot the antenna locations
% scatter(x_locations, y_locations, 'r', 'filled');
% 
% % Optionally, add labels to the antenna locations
% for i = 1:length(x_locations)
%     text(x_locations(i), y_locations(i), sprintf('A%d', i), 'Color', 'white', 'FontSize', 12, 'HorizontalAlignment', 'center');
% end
% 
% % Release the hold on the current plot
% hold off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THRESHOLD %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set the threshold value
threshold = 0.075;
% Create a copy of the original image slice
im_slice_1 = im_slice;

% Perform threshold segmentation
for i = 1:size(im_slice_1, 1)
    for j = 1:size(im_slice_1, 2)
        if im_slice_1(i, j) < threshold
            im_slice_1(i, j) = 0; % Set values below the threshold to 0
        end
    end
end

% Display the segmented image
imagesc(slice_axes_x, slice_axes_y, im_slice_1');
xlabel('X-axis');
ylabel('Y-axis');

% Calculate the number of non-zero elements after thresholding
sum_t = 0;
for i = 1:size(im_slice_1, 1)
    for j = 1:size(im_slice_1, 2)
        if im_slice_1(i, j) > 0 && ~isnan(im_slice_1(i, j))
            sum_t = sum_t + 1;
        end
    end
end
disp(sum_t);

% Calculate the number of non-zero elements in the original image
sum_non_0 = 0;
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if im_slice(i, j) > 0 && ~isnan(im_slice(i, j))
            sum_non_0 = sum_non_0 + 1;
        end
    end
end
disp(sum_non_0);

% Calculate and display the percentage of non-zero elements retained
percent = sum_t / sum_non_0 * 100;
disp([num2str(percent) '%']);

% Initialize counters for different threshold ranges
sum_1 = 0;
sum_2 = 0;
sum_3 = 0;

% Count the number of elements above 0.05
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if im_slice(i, j) > 0.05 && ~isnan(im_slice_1(i, j))
            sum_1 = sum_1 + 1;
        end
    end
end

% Count the number of elements between 0.025 and 0.05
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if im_slice(i, j) > 0.025 && im_slice(i, j) < 0.05 && ~isnan(im_slice_1(i, j))
            sum_2 = sum_2 + 1;
        end
    end
end

% Count the number of elements below 0.025
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if im_slice(i, j) < 0.025 && ~isnan(im_slice_1(i, j))
            sum_3 = sum_3 + 1;
        end
    end
end

% Display the counts for each range
disp(sum_1);
disp(sum_2);
disp(sum_3);

% Calculate the total sum of all elements in the original image
sum_in = 0;
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if ~isnan(im_slice(i, j))
            sum_in = sum_in + im_slice(i, j);
        end
    end
end
disp(sum_in);
disp(z_slice_value);

% Calculate and display the peak-to-peak value in dB
db = 10 * log10(sum_in);
X = ['db=', num2str(db)];
disp(X);

%%%%%%%%%%%%%%% CODE #1 %%%%%%%%%%%%%%%
% % Initialize sums and counts for signal and noise
% sum_s = 0;
% sum_n = 0;
% count_s = 0;
% count_n = 0;
% 
% % Loop through each element of the image slice
% for i = 1:size(im_slice, 1)
%     for j = 1:size(im_slice, 2)
%         if im_slice(i, j) > threshold && ~isnan(im_slice(i, j))
%             sum_s = sum_s + im_slice(i, j); % Sum of signal values
%             count_s = count_s + 1; % Count of signal values
%         elseif im_slice(i, j) <= threshold && ~isnan(im_slice(i, j))
%             sum_n = sum_n + im_slice(i, j); % Sum of noise values
%             count_n = count_n + 1; % Count of noise values
%         end    
%     end
% end
% 
% % Calculate the average signal and noise values
% if count_s > 0
%     average_s = sum_s / count_s;
% else
%     average_s = 0;
% end
% 
% if count_n > 0
%     average_n = sum_n / count_n;
% else
%     average_n = 0;
% end
% 
% % Initialize variances for signal and noise
% variance_s = 0;
% variance_n = 0;
% 
% % Calculate variance for signal and noise
% for i = 1:size(im_slice, 1)
%     for j = 1:size(im_slice, 2)
%         if im_slice(i, j) > threshold && ~isnan(im_slice(i, j))
%             variance_s = variance_s + (im_slice(i, j) - average_s)^2; % Signal variance
%         elseif im_slice(i, j) <= threshold && ~isnan(im_slice(i, j))
%             variance_n = variance_n + (im_slice(i, j) - average_n)^2; % Noise variance
%         end    
%     end
% end
% 
% % Normalize the variances - the variances for signal and noise are normalized 
% % by dividing by count_s - 1 for signal and count_n - 1 for noise
% if count_s > 1
%     variance_s = variance_s / (count_s - 1);
% else
%     variance_s = 0;
% end
% 
% if count_n > 1
%     variance_n = variance_n / (count_n - 1);
% else
%     variance_n = 0;
% end

% % Calculate SNR: signal power over noise power
% if variance_n > 0
%     snr = 10 * log10(variance_s / variance_n); % SNR in dB
% else
%     snr = Inf; % If noise variance is 0, SNR is infinite
% end
% 
% % Display SNR
% X = ['snr=', num2str(snr)];
% disp(X);

%%%%%%%%%%%%%%%%%%%%%%%%%% CODE #2 %%%%%%%%%%%%% SQUARING PIXEL %%%%%%%%%%%
% Initialize variables
sum_s = 0;
sum_n = 0;
count_s = 0;
count_n = 0;

% Iterate over the image slice
for i = 1:size(im_slice, 1)
    for j = 1:size(im_slice, 2)
        if im_slice(i, j) > threshold && ~isnan(im_slice(i, j))
            sum_s = sum_s + im_slice(i, j)^2; % Sum of squared signal pixel intensities
            count_s = count_s + 1; % Count of signal pixels
        elseif im_slice(i, j) <= threshold && ~isnan(im_slice(i, j))
            sum_n = sum_n + im_slice(i, j)^2; % Sum of squared noise pixel intensities
            count_n = count_n + 1; % Count of noise pixels
        end    
    end
end

% Calculate the mean square values (power) of signal and noise
mean_square_s = sum_s / count_s;
mean_square_n = sum_n / count_n;

% Calculate and display SNR in dB
snr = 10 * log10(mean_square_s / mean_square_n);
X = ['SNR = ', num2str(snr), ' dB'];
disp(X);

