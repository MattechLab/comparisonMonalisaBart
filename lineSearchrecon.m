function lineSearchrecon(prefix, regtype, regvals)
    % LINESEARCHRECON - Find the best regularization parameters for L1 and L2
    % regularized reconstructions. Save the values and plot the curves.
    %
    % inputs: prefix (can be ['eye','phantom','cardiac'] identifies which
    % image to use), regtype (can be ['l1','l2'] identifies the type of
    % regularization in the reconstruction), regvals (list of regualarization values to test)
    %This function will:
    % 1. Generate the simulated data
    % 2. Run the reconstruction for many different regularization values,
    % all of them are contained in regvals
    % 3. Normalize the reconstruction results on the ROI, so that it
    % matches the groud truth
    % 4. Plot l2distances and ssim, to identify the best regval.


    %% accquisition and recon parameters
    
    N       = 256; % Size of the image (e.g., 256x256)
    N_u     = [256, 256];
    n_u     = N_u;
    nLines  = 30;
    
    % Select the ground truth image and ellipse based on prefix
    % Select the coil sense accordingly
    % Create an ellipse on the ROI 
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
    
    dK_u    = 1./FoV;
    % 1. SIMULATED TRAJECTORY t_tot
    t_tot = bmTraj_fullRadial2_lineAssym2(N,nLines, dK_u(1));
    
    % 2. SIMULATED DATA y
    
    % Total number of points
    numPoints = size(t_tot, 2);
    
    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
    
    % 3. RUN MATHILDA: OPTIMALLY SELECT A NOT SOO NICE MATHILDA THAT IMPROVES
    
    ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
    
    % GRIDDED RECON
    x0 = bmMathilda(y,t_tot,ve,C, N_u, n_u, dK_u);
    
    % SAVE GRIDDED RECON


    % Gridding matrices for iterative recons
    [Gu, Gut] = bmTraj2SparseMat(t_tot, ve, N_u, dK_u);

    %% PARAMS FOR THE REG reconstruction with monalisa
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
            frameSize = n_u;
            rho = 10*delta; % How to select rho?? => rule of thumb 10*delta also to bart
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
        image_mean_roi = mean(abs(image(roi_mask))); % gt mean value
        x_mean_roi = mean(abs(x(roi_mask))); % reconstruction mean value
        
        % Apply normalization if ROI mean of reconstruction is nonzero
        x = x * (image_mean_roi / x_mean_roi);
        pred = x;
        % Allign the ground truth to the prediction to compute the metrics
        % New: align ground truth magnitude to prediction magnitude
        [pred_aligned_gt, a, b] = align_gt_to_recon_magnitude(image, x);

        % Compute L2 distance between ground truth image and predicted image
        %diff = image(:) - pred(:);
        pred_mag = abs(pred);
        aligned_gt_mag = abs(pred_aligned_gt);

        diff_mags = (aligned_gt_mag(:) - pred_mag(:));
        l2distance = norm(diff_mags);  % equivalent to sqrt(sum(diff.^2))

        %diff = pred_aligned_gt(:) - pred(:);
        %l2distance = sqrt(real(diff' * diff)); % we compute the l2 distance of the magnitude because of the alignment of the ground truth

        % Append the L2 distance and regularization parameter to the arrays
        l2_distances = [l2_distances, l2distance];
    
        % Compute SSIM of the magnitude of the prediction w.r.t. the
        % magnitude of the gt
        
        % Compute robust dynamic range from percentiles
        lo = prctile(aligned_gt_mag(:), 0.5);
        hi = prctile(aligned_gt_mag(:), 99.5);
        L = hi - lo;
        % SSIM with explicit dynamic range
        ssim_value = ssim(double(pred_mag), double(aligned_gt_mag), 'DynamicRange', L);
        %image_mag = abs(image);
        %ssim_value = ssim(pred_mag, image_mag);

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

