% generate the simulated rawdata for the comparison.
% FOR BART RECONSTUCTION WE NEED:
% 0. the ground truth files 
% 1. A trajectory file
% 2. A volume element file
% 3. A raw data file
% 4. A coil sensitivity file
% Each file should be saved in cfl file format, we can use the
% helperfunction writecfl. NB: The dimentions of the array have a meaning.
% Refear to Bart documentation to further understand bart conventions

clear;
close all;
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

%% SL-Phantom
bmImage(shepp_logan_phantom);

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
bmImage(finalbrainimage)

%% Cardiac image
cardiacimage = abs(h{1});
bmImage(cardiacimage);

cardiacimage2 = abs(single(imresize(x{1}, [512, 512])));
bmImage(cardiacimage2);

cardiacimage3 = load('./images/h_3.mat');
cardiacimage3 = abs(single(imresize(cardiacimage3.h, [512, 512])));

%% SAVE IN BART FORMAT: naming convention is final + prefix + image
writecfl('bart_data/finalphantomimage', shepp_logan_phantom);
writecfl('bart_data/finalbrainimage', finalbrainimage);
writecfl('bart_data/finalcardiacimage', cardiacimage);
writecfl('bart_data/finalcardiac2image', cardiacimage2);
writecfl('bart_data/finalcardiac3image', cardiacimage3);
C = bmImResize(C, [64, 64], N_u);
% Convert and save Coil Sensitivities (C)
C_reshaped = reshape(C, [512, 512, 1, 30]);
% Convert and save Coil Sensitivities (C)
writecfl('bart_data/C_bart', C_reshaped);

% Define image-prefix pairs
image_prefix_pairs = {...
shepp_logan_phantom, 'phantom'; ...
finalbrainimage, 'brain'; ...
cardiacimage, 'cardiac'; ...
cardiacimage2, 'cardiac2'; ...
cardiacimage3, 'cardiac3'};

% 1. SIMULATED TRAJECTORY t_tot
nLines = 30;
t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));

% 2. SIMULATED DATA y
numPoints = size(t_tot, 2);



% 3. MATHILDA PROCESSING
ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
rescaler = FoV(1);
traj_rescaled = t_tot * rescaler;
traj_reshaped = reshape(traj_rescaled, [2, N, nLines]);
traj_bart = cat(1, traj_reshaped, zeros(1, N, nLines));

% Validate reshaping
idx = randi([1, N*nLines]);
assert(isequal(traj_rescaled(:,idx), traj_reshaped(:, mod(idx-1, 512) + 1, floor((idx-1) / 512) + 1)))

% Write trajectory data
writecfl('bart_data/traj_bart', traj_bart);

ve_max = 3/prod(FoV(:));
ve_srt_cap = sqrt(min(ve, ve_max));
ve_srt_cap = ve_srt_cap*prod(FoV(:));
% Convert and save Volume Elements (weights)
weights_bart = reshape(ve_srt_cap, [1, N, nLines]);

assert(isequal(weights_bart(1,3,5),ve_srt_cap(1,(5-1)*N + 3)))
writecfl('bart_data/weights_bart', weights_bart);

% Loop over each image-prefix pair
for i = 1:size(image_prefix_pairs, 1)
    image = image_prefix_pairs{i, 1};
    prefix = image_prefix_pairs{i, 2};

    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
    % Convert and save k-space data
    kspace_reshaped = reshape(y, [N, nLines, 30]);
    kspace_bart = reshape(kspace_reshaped, [1, N, nLines, 30, 1]);
    
    % Validate reshaping
    idx2 = randi([1, 30]);
    assert(isequal(y(idx,idx2), kspace_bart(:, mod(idx-1, 512) + 1, floor((idx-1) / 512) + 1, idx2, :)))
    
    % Write k-space data
    writecfl(['bart_data/kspace_bart', prefix], kspace_bart);

    % Run gridded recon (Mathilda), and save the cfl files
    % GRIDDED RECON:
    griddedrecon = bmMathilda(y,t_tot,ve_srt_cap,C, N_u, n_u, dK_u);
    writecfl(['reconstructions/griddedrecons/gridded_recon_final_', prefix], griddedrecon);
end


    







