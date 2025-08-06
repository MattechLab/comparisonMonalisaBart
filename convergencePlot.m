function [] = convergencePlot(prefix,regtype)
%convergencePlot Given a wintessInfo object and the associated ground truth
%image generate a convergence plot for the iterative reconstruction.
%   Given a wintessInfo object and the associated ground truth
%    image generate a convergence plot for the iterative reconstruction.
witnessInfoPath = ['./witnessInfos/wit_monalisa_final_recon_',prefix,'_',regtype '.mat'];
N = 256;
switch prefix
        case 'phantom'
            load('./images/C_simu.mat'); % Load the simulated coil sensitivity for the phantom
            C = bmImResize(C, [96, 96], [N,N]); % fix C to same dimentions
            image = single(phantom('Modified Shepp-Logan', N));
            ellipse = [129, 127, 88, 117.5, 0]; % Ellipse on the phantom
            FoV     = [240, 240]; % Random FoV for simulation with phantom
        
        case 'eye'
            load('./images/image_eye_yiwei.mat'); % Load the brain/eye image from yiwei and the C
            x = bmImResize(x, [480, 480], [N,N]); % downsize x to 256
            image = single(x); 
            C = bmImResize(C, [480, 480], [N,N]); % downsize C to 256
            ellipse = [127, 85, 85, 110, 0]; % Ellipse on the brain
            FoV     = [240, 240]; % double FoV yiwei

        case 'cardiac'
            load('./images/image_cardiac_bern.mat')
            image = single(x);
            ellipse = [143, 132.5, 64, 45, 68]; % Check if it's correct
            FoV     = [300, 300]; % double FoV cardiac acquisition
            
        otherwise
            error('Invalid prefix provided. The image does not exist');% Ellipse on the hearth
end

% Initialize empty arrays to store results
l2_distances = [];  % For L2 distance
ssim_values = [];   % For SSIM

% Load the witnessinfo images
load(witnessInfoPath);
witnessImgs = witnessInfo.witness_im;
% Compute the timeseries values of scores for l2 distance and SSIM
for i = 1:size(witnessImgs,3)
    x = witnessImgs(:,:,i);
    % Normalize reconstruction over ROI using magnitude (abs)
    [xx, yy] = meshgrid(1:N, 1:N);
    center_x = ellipse(1); center_y = ellipse(2);
    major_axis = ellipse(3); minor_axis = ellipse(4);
    theta = deg2rad(ellipse(5));
    x_rot = cos(theta) * (xx - center_x) + sin(theta) * (yy - center_y);
    y_rot = -sin(theta) * (xx - center_x) + cos(theta) * (yy - center_y);
    roi_mask = ((x_rot / major_axis).^2 + (y_rot / minor_axis).^2) <= 1;
    
    % Compute mean of the magnitude (abs) within the ROI
    image_mean_roi = mean(abs(image(roi_mask))); % gt mean value
    x_mean_roi = mean(abs(x(roi_mask))); % reconstruction mean value
    
    % Apply normalization if ROI mean of reconstruction is nonzero
    x = x * (image_mean_roi / x_mean_roi);
    pred = x;
    % Compute L2 distance between ground truth image and predicted image
    diff = image(:) - pred(:);
    l2distance = sqrt(real(diff' * diff)); % handle complex values

    % Append the L2 distance and regularization parameter to the arrays
    l2_distances = [l2_distances, l2distance];

    % Compute SSIM of the magnitude of the prediction w.r.t. the
    % magnitude of the gt
    pred_mag = abs(pred);
    image_mag = abs(image);
    ssim_value = ssim(pred_mag, image_mag);
    ssim_values = [ssim_values, ssim_value];
end

% Rescale for better plotting?
ssim_values = (ssim_values-min(ssim_values))/(max(ssim_values-min(ssim_values)));
l2_distances = (l2_distances-min(l2_distances))/(max(l2_distances)-min(l2_distances));
% Make the plot
plot(1:size(witnessImgs,3),ssim_values,1:size(witnessImgs,3),l2_distances)
title(sprintf('Convergence of SSIM and L2 Distance for %s (%s Regularization)', prefix, upper(regtype)));
xlabel('Iteration number')
legend('structural similarity','l2 distance (min-max rescaled)')
end