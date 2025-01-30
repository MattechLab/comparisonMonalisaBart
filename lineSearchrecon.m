function lineSearchrecon(prefix, regtype, regvals)
    % LINESEARCHRECON - Find the best regularization parameters for L1 and L2
    % regularized reconstructions. Save the values and plot the curves.
    load('./images/x_cartesian.mat')
    % add monalisa to path
    addpath(genpath('/Users/mauroleidi/Desktop/monalisa/src')); 
    
    % read the ground truth image (NB: this has to be replaced afterwards)
    cinefilepaht = './images/radial_2D_CINE_bastien.mat';
    load(cinefilepaht);
    
    % Define the size of the phantom image
    N = 512; % Size of the image (e.g., 256x256)
    % Generate the 2D Shepp-Logan phantom
    shepp_logan_phantom = phantom('Shepp-Logan', N);
    shepp_logan_phantom = single(min(10*shepp_logan_phantom,1));
    
    %% 2D brain image
    
    % Load the PNG image
    brainimage = imread('./images/brain_image.png');
    % Convert to grayscale
    gray_brainimage = rgb2gray(brainimage);
    % Crop the center (367x367) from the original image
    shift = 11;
    crop_size = 367;
    center_x = round(size(gray_brainimage, 2) / 2) - shift; % X center of the image
    center_y = round(size(gray_brainimage, 1) / 2); % Y center of the image
    % Crop the image around the center
    cropped_brainimage = gray_brainimage(center_y - floor(crop_size/2) + 1:center_y + floor(crop_size/2), ...
                                          center_x - floor(crop_size/2) + 1:center_x + floor(crop_size/2));
    % Resize the cropped image to 512x512
    finalbrainimage = single(imresize(cropped_brainimage, [512, 512]));
    
    
    % Cardiac images
    cardiacimage = abs(h{1});    
    cardiacimage2 = abs(single(imresize(x{1}, [512, 512])));
   
    % Select the image and ellipse based on prefix
    switch prefix
        case 'phantom'
            image = shepp_logan_phantom;
            ellipse = [258, 254, 176, 235, 0]; 
        case 'brain'
            image = finalbrainimage;
            ellipse = [271, 250, 202, 220, 0];
        case 'cardiac'
            image = cardiacimage;
            ellipse = [250, 234, 107, 95, 75];
        case 'cardiac2'
            image = cardiacimage2;
            ellipse = [262, 217, 100, 132, 145];
        otherwise
            error('Invalid prefix provided.');
    end

    % 1. SIMULATED TRAJECTORY t_tot
    nLines = 30;
    t_tot = bmTraj_fullRadial2_lineAssym2(N,nLines, dK_u(1));
    
    % 2. SIMULATED DATA y
    
    % Total number of points
    numPoints = size(t_tot, 2);
    
    C = bmImResize(C, [64, 64], N_u);
    
    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
    
    % 3. RUN MATHILDA: OPTIMALLY SELECT A NOT SOO NICE MATHILDA THAT IMPROVES
    
    ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
    
    % GRIDDED RECON: Warm start
    x0 = bmMathilda(y,t_tot,ve,C, N_u, n_u, dK_u);
    
    [Gu, Gut] = bmTraj2SparseMat(t_tot, ve, N_u, dK_u);
    %% PARAMS FOR THE L2 REG reconstruction with monalisa
    nIter = 35;
    nCGD = 4;
    witness_ind   = 1:5:nIter; % indices to save the witness
    witness_label = 'label'; % label to save the file
    regul_mode = 'normal';
    ve_max = 3/prod(FoV(:));

    
    % Initialize empty arrays to store results
    l2_distances = [];  % For L2 distance
    ssim_values = [];   % For SSIM
    
    % Loop over the regularization values
    for i = 1:length(regvals)
        % Set the current regularization value
        delta = regvals(i);
        
        if strcmp(regtype, 'l2')
            x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, delta, regul_mode, nCGD, ve_max, ...
                        nIter, bmWitnessInfo(witness_label, witness_ind));
        elseif strcmp(regtype, 'l1')
            % How to select rho?? => rule of thumb 10*delta also to bart
            frameSize = n_u;
            rho = 10*delta; % this should be some what adapted too.
            x = bmSteva(x0, [], [], y, ve, C, Gu, Gut, frameSize, delta, rho, nCGD, ve_max, ...
                        nIter, bmWitnessInfo(witness_label, witness_ind));
        end
        
        % Normalize reconstruction over ROI using magnitude (abs)
        [xx, yy] = meshgrid(1:N, 1:N);
        center_x = ellipse(1); center_y = ellipse(2);
        major_axis = ellipse(3); minor_axis = ellipse(4);
        theta = deg2rad(ellipse(5));
        x_rot = cos(theta) * (xx - center_x) + sin(theta) * (yy - center_y);
        y_rot = -sin(theta) * (xx - center_x) + cos(theta) * (yy - center_y);
        roi_mask = ((x_rot / major_axis).^2 + (y_rot / minor_axis).^2) <= 1;
        
        % Compute mean of the magnitude (abs) within the ROI
        image_mean_roi = mean(abs(image(roi_mask)));
        x_mean_roi = mean(abs(x(roi_mask)));
        
        % Apply normalization if ROI mean of reconstruction is nonzero
        x = x * (image_mean_roi / x_mean_roi);
    
        % Get the model prediction (assuming x is the predicted image)
        pred = abs(x);  % Ensure absolute value if it's a complex image
    
        % Compute L2 distance between ground truth image and predicted image
        l2distance = sqrt(sum((image(:) - pred(:)).^2));
        % Append the L2 distance and regularization parameter to the arrays
        l2_distances = [l2_distances, l2distance];
    
        % Compute SSIM
        ssim_value = ssim(pred, image);
        ssim_values = [ssim_values, ssim_value];
    
        
    end
    
    [maxValue, index] = max(ssim_values);
    bestregval = regvals(index);
    
    % Plot all metrics on the same figure
    figure;
    yyaxis left; % Left y-axis for L2 Distance
    plot(regvals, l2_distances, '-o', 'LineWidth', 1.5, 'DisplayName', 'L2 Distance');
    ylabel('L2 Distance');
    set(gca, 'YColor', 'b');  % Set color for left y-axis
    
    yyaxis right; % Right y-axis for SSIM and PSNR
    hold on;
    plot(regvals, ssim_values, '-s', 'LineWidth', 1.5, 'DisplayName', 'SSIM');
    ylabel('SSIM');
    set(gca, 'YColor', 'r');  % Set color for right y-axis
    
    % Common x-axis
    xlabel('Regularization Parameter');
    set(gca, 'XScale', 'log');  % Log scale for regularization parameter
    % Dynamically set title to include regularization type
    title(['Evaluation Metrics Across Regularization Parameters (', upper(regtype), ')', prefix, 'best regval: ',num2str(bestregval)]);
    grid on;
    
    % Add legend
    legend('show', 'Location', 'best');

    % SAVE RESULTING IMAGES IN SOME FOLDER TOGHETHER WITH ALL THE SSIM VALS
    % AND REGVALS AND TYPE OF REG AND PREFIX
    %% Save Results
    %results_path = 'results_monalisa/';
    %if ~exist(results_path, 'dir')
    %    mkdir(results_path);
    %end
    %save(fullfile(results_path, ['results_', prefix, '_', regtype, '.mat']), 'regvals', 'l2_distances', 'ssim_values', 'bestregval');
end

