function x = runFinalReconstruction(prefix, regtype, bestregval, saveFlag)
    % RUNFINALRECONSTRUCTION - Run and optionally save final reconstruction using optimal regularization.
    %
    % inputs: 
    %   prefix    - ['eye','phantom','cardiac'] identifies which image to use
    %   regtype   - ['l1','l2'] identifies the type of regularization
    %   bestregval- optimal regularization value from gridsearch
    %   saveFlag  - (optional) boolean, whether to save results (default = true)
    %
    %This function will:
    % 1. Generate the simulated data
    % 2. Run the reconstruction and normalize the result on a ROI
    % 3. Store the witness indicator to make the convergence plots
    % 4. Optionally save the reconstructions in ./reconstructions/monalisa/

    if nargin < 4
        saveFlag = true; % default
    end

    %% acquisition and recon parameters
    N       = 256; % Size of the image (e.g., 256x256)
    N_u     = [256, 256];
    n_u     = N_u;
    nLines  = 30;
    
    % Select the ground truth image and ellipse based on prefix
    switch prefix
        case 'phantom'
            load('./images/C_simu.mat'); 
            C = bmImResize(C, [96, 96], [N,N]);
            image = phantom('Modified Shepp-Logan', N);
            ellipse = [129, 127, 88, 117.5, 0];
            FoV     = [240, 240];
        
        case 'eye'
            load('./images/image_eye_yiwei.mat');
            x = bmImResize(x, [480, 480], [N,N]);
            image = x; 
            C = bmImResize(C, [480, 480], [N,N]);
            ellipse = [127, 85, 85, 110, 0];
            FoV     = [240, 240];

        case 'cardiac'
            load('./images/image_cardiac_bern.mat')
            image = x;
            ellipse = [143, 132.5, 64, 45, 68];
            FoV     = [300, 300];
            
        otherwise
            error('Invalid prefix provided.');
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

    % Witness info
    witness_ind   = 1:5:nIter;
    witness_label = ['monalisa_final_recon_',prefix,'_',regtype];
    witness_info = bmWitnessInfo(witness_label, witness_ind);
    witness_info.save_witnessIm_flag = saveFlag;

    % Run reconstruction
    if strcmp(regtype, 'l2')
        x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, bestregval, regul_mode, nCGD, ve_max, nIter, witness_info);
    elseif strcmp(regtype, 'l1')
        rho = 10 * bestregval;
        x = bmSteva(x0, [], [], y, ve, C, Gu, Gut, n_u, bestregval, rho, nCGD, ve_max, nIter, witness_info);
    else
        error('Invalid regularization type.');
    end

    % Normalize over ROI
    [xx, yy] = meshgrid(1:N, 1:N);
    center_x = ellipse(1); center_y = ellipse(2);
    major_axis = ellipse(3); minor_axis = ellipse(4);
    theta = deg2rad(ellipse(5));
    x_rot = cos(theta) * (xx - center_x) + sin(theta) * (yy - center_y);
    y_rot = -sin(theta) * (xx - center_x) + cos(theta) * (yy - center_y);
    roi_mask = ((x_rot / major_axis).^2 + (y_rot / minor_axis).^2) <= 1;
    
    image_mean_roi = mean(abs(image(roi_mask)));
    x_mean_roi = mean(abs(x(roi_mask)));
    x = x * (image_mean_roi / x_mean_roi);

    % Save results (only if saveFlag is true)
    if saveFlag
        results_folder = './reconstructions/monalisa/';
        if ~exist(results_folder, 'dir')
            mkdir(results_folder);
        end
        save_path = fullfile(results_folder, sprintf('%s_%s.mat', prefix, regtype));
        save(save_path, 'x', 'bestregval', 'prefix', 'regtype');
        fprintf('Reconstruction saved to %s\n', save_path);
    else
        fprintf('Reconstruction not saved (saveFlag=false).\n');
    end
end
