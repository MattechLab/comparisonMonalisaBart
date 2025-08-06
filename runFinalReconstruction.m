function runFinalReconstruction(prefix, regtype, bestregval)
    % RUNFINALRECONSTRUCTION - Run and save final reconstruction using optimal regularization.
    % inputs: prefix (can be ['eye','phantom','cardiac'] identifies which
    % image to use), regtype (can be ['l1','l2'] identifies the type of
    % regularization in the reconstruction), bestregval (the optimal
    % regularization value estimated using a gridsearch) 
    %
    %This function will:
    % 1. Generate the simulated data
    % 2. Run the reconstruction and normalize the result on a ROI
    % 3. Store the witness indicator to make the convergence plots
    % 4. Save the reconstuctions in the path ./reconstructions/monalisa/prefix_regtype.mat

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
            image = phantom('Modified Shepp-Logan', N);
            ellipse = [129, 127, 88, 117.5, 0]; % Ellipse on the phantom
            FoV     = [240, 240]; % Random FoV for simulation with phantom
        
        case 'eye'
            load('./images/image_eye_yiwei.mat'); % Load the brain/eye image from yiwei and the C
            x = bmImResize(x, [480, 480], [N,N]); % downsize x to 256
            image = x; 
            C = bmImResize(C, [480, 480], [N,N]); % downsize C to 256
            ellipse = [127, 85, 85, 110, 0]; % Ellipse on the brain
            FoV     = [240, 240]; % double FoV yiwei

        case 'cardiac'
            load('./images/image_cardiac_bern.mat')
            image = x;
            ellipse = [143, 132.5, 64, 45, 68]; % Check if it's correct
            FoV     = [300, 300]; % double FoV cardiac acquisition
            
        otherwise
            error('Invalid prefix provided. The image does not exist');% Ellipse on the hearth
    end
    
    dK_u    = 1./FoV;
    
    t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
    ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
    x0 = bmMathilda(y, t_tot, ve, C, N_u, n_u, dK_u);
    [Gu, Gut] = bmTraj2SparseMat(t_tot, ve, N_u, dK_u);
    
    % Reconstruction parameters
    nIter = 160;
    nCGD = 4;
    regul_mode = 'normal';
    ve_max = 3/prod(FoV(:));

    % Setting up witness indicator to check convergence
    witness_ind   = 1:5:nIter; % indices to save the witness
    witness_label = ['monalisa_final_recon_',prefix,'_',regtype]; % label to save the file
    witness_info = bmWitnessInfo(witness_label, witness_ind);
    witness_info.save_witnessIm_flag = true;
    % Run final reconstruction
    if strcmp(regtype, 'l2')
        x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, bestregval, regul_mode, nCGD, ve_max, nIter, witness_info);
    elseif strcmp(regtype, 'l1')
        rho = 10 * bestregval;
        x = bmSteva(x0, [], [], y, ve, C, Gu, Gut, n_u, bestregval, rho, nCGD, ve_max, nIter, witness_info);
    else
        error('Invalid regularization type.');
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

    % Save results
    results_folder = './reconstructions/monalisa/';
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
    save_path = fullfile(results_folder, sprintf('%s_%s.mat', prefix, regtype));
    save(save_path, 'x', 'bestregval', 'prefix', 'regtype');
    fprintf('Reconstruction saved to %s\n', save_path);
end
