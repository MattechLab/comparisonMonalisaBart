% generate the simulated rawdata for the comparison.
% FOR BART RECONSTUCTION WE NEED:
% 0. the ground truth files 
% 1. A trajectory file (common for all the recon)
% 2. A volume element filem (common for all the recon)
% 3. A raw data file (simulated from the ground truth)
% 4. A coil sensitivity file
% Each file should be saved in cfl file format, we can use the
% helperfunction writecfl. NB: The dimentions of the array have a meaning.
% Refear to Bart documentation to further understand bart conventions and
% why we reformat arrays before writing them to cfl

% you can also check our discussions here https://github.com/mrirecon/bart/issues/345
% BIG THANKS TO THE BART TEAM!! 

clear;
close all;

% add monalisa to path
addpath(genpath('/Users/mauroleidi/Desktop/monalisa/src')); 

% Define prefixes
image_prefixes = {'phantom';'eye';'cardiac'};

% Acquisition and simulation parameters
N = 256; % Size of the image (e.g., 256x256)
nLines = 30;
N_u     = [256, 256];
n_u     = N_u;

% Loop over each prefixes and generate the Bart files
for i = 1:size(image_prefixes)
    prefix = image_prefixes{i};
    disp(['Currently saving bart data for: ',prefix])
    % Switch to select image and C based on the prefix
    switch prefix
        case 'phantom'
            load('./images/C_simu.mat'); % Load the simulated coil sensitivity for the phantom
            C = bmImResize(C, [96, 96], [N,N]); % fix C to same dimentions
            image = single(phantom('Modified Shepp-Logan', N));
            FoV     = [240, 240]; % Random FoV for simulation with phantom
        
        case 'eye'
            load('./images/image_eye_yiwei.mat'); % Load the brain/eye image from yiwei and the C
            x = bmImResize(x, [480, 480], [N,N]); % downsize x to 256
            image = single(x); 
            C = bmImResize(C, [480, 480], [N,N]); % downsize C to 256
            FoV     = [240, 240]; % double FoV yiwei

        case 'cardiac'
            load('./images/image_cardiac_bern.mat')
            image = single(x);
            FoV     = [300, 300]; % double FoV cardiac acquisition
    end

    %% SAVE IMAGE IN BART FORMAT: naming convention is final + prefix + image
    writecfl(['bart_data/image_bart_', prefix], image);
    
    nCh = size(C,3);
    % Convert and save Coil Sensitivities (C)
    C_reshaped = reshape(C, [N, N, 1, nCh]);
    %% SAVE C IN BART FORMAT: naming convention is C_bart + prefix
    writecfl(['bart_data/C_bart_',prefix], C_reshaped);

    dK_u    = 1./FoV;

    %% let's compute the trajectory and volume elements 
    % 1. SIMULATED TRAJECTORY t_tot
    t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
    
    % 2. SIMULATED DATA y
    numPoints = size(t_tot, 2);
    
    % 3. MATHILDA PROCESSING
    rescaler = FoV(1);
    traj_rescaled = t_tot * rescaler;
    traj_reshaped = reshape(traj_rescaled, [2, N, nLines]);
    traj_bart = cat(1, traj_reshaped, zeros(1, N, nLines));
    
    % Validate reshaping
    idx = randi([1, N*nLines]);
    assert(isequal(traj_rescaled(:,idx), traj_reshaped(:, mod(idx-1, N) + 1, floor((idx-1) / N) + 1)))
    
    %% SAVE TRAJ IMAGE IN BART FORMAT
    writecfl(['bart_data/traj_bart_',prefix], traj_bart);
    
    % Compute volume elements
    ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
    ve_max = 3/prod(FoV(:));
    ve = min(ve, ve_max);
    % Restructure following BART conventions
    ve_srt_cap = sqrt(ve);
    % This is made to follow the BART convention about density compensation.
    ve_srt_cap = ve_srt_cap*prod(FoV(:)); 
    % Convert and save Volume Elements (weights)
    weights_bart = reshape(ve_srt_cap, [1, N, nLines]);
    % Validate reshaping
    assert(isequal(weights_bart(1,3,5),ve_srt_cap(1,(5-1)*N + 3)))
    
    writecfl(['bart_data/weights_bart_',prefix], weights_bart);

    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);

    % Convert and save k-space data
    kspace_reshaped = reshape(y, [N, nLines, nCh]);
    kspace_bart = reshape(kspace_reshaped, [1, N, nLines, nCh, 1]);
    
    % Validate reshaping
    idx2 = randi([1, nCh]);
    assert(isequal(y(idx,idx2), kspace_bart(:, mod(idx-1, N) + 1, floor((idx-1) / N) + 1, idx2, :)))
    
    % Write k-space data
    writecfl(['bart_data/kspace_bart_', prefix], kspace_bart);

    % Run gridded recon (Mathilda), and save the cfl files
    % GRIDDED RECON FOR WARM START:
    griddedrecon = bmMathilda(y,t_tot,ve,C, N_u, n_u, dK_u);
    writecfl(['reconstructions/griddedrecons/gridded_recon_', prefix], griddedrecon);
end
