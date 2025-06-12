function runFinalReconstruction(prefix, regtype, bestregval)
    % RUNFINALRECONSTRUCTION - Run and save final reconstruction using optimal regularization.
    
    % Load necessary data
    load('./images/x_cartesian.mat');
    cinefilepaht = './images/radial_2D_CINE_bastien.mat';
    load(cinefilepaht);
    addpath(genpath('/Users/mauroleidi/Desktop/monalisa/src'));
    
    % Define image selection based on prefix
    N = 512;
    shepp_logan_phantom = phantom('Shepp-Logan', N);
    shepp_logan_phantom = single(min(10*shepp_logan_phantom,1));
    brainimage = imread('./images/brain_image.png');
    gray_brainimage = rgb2gray(brainimage);
    crop_size = 367;
    shift = 11;
    center_x = round(size(gray_brainimage, 2) / 2) - shift;
    center_y = round(size(gray_brainimage, 1) / 2);
    cropped_brainimage = gray_brainimage(center_y - floor(crop_size/2) + 1:center_y + floor(crop_size/2), ...
                                          center_x - floor(crop_size/2) + 1:center_x + floor(crop_size/2));
    finalbrainimage = single(imresize(cropped_brainimage, [512, 512]));
    cardiacimage = abs(h{1});
    cardiacimage2 = abs(single(imresize(x{1}, [512, 512])));
    
    cardiacimagenew = load('/Users/mauroleidi/Desktop/comparisonMonalisaBartNew/comparisonMonalisaBart/images/h_3.mat');
    cardiacimagenew = abs(single(imresize(cardiacimagenew.h, [512, 512])));

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
        case 'cardiacnew'
            image = cardiacimagenew;
            ellipse = [273, 268, 128, 128, 0];
        otherwise
            error('Invalid prefix provided.');
    end
    
    % Generate trajectory and data
    nLines = 30;
    t_tot = bmTraj_fullRadial2_lineAssym2(N, nLines, dK_u(1));
    C = bmImResize(C, [64, 64], N_u);
    y = bmSimulateMriData(image, C, t_tot, N_u, n_u, dK_u);
    ve = bmVolumeElement(t_tot, 'voronoi_full_radial2');
    x0 = bmMathilda(y, t_tot, ve, C, N_u, n_u, dK_u);
    [Gu, Gut] = bmTraj2SparseMat(t_tot, ve, N_u, dK_u);
    
    % Reconstruction parameters
    nIter = 35;
    nCGD = 4;
    regul_mode = 'normal';
    ve_max = 3/prod(FoV(:));
    witness_ind   = {}; % indices to save the witness
    witness_label = 'label'; % label to save the file

    % Run final reconstruction
    if strcmp(regtype, 'l2')
        x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, bestregval, regul_mode, nCGD, ve_max, nIter, bmWitnessInfo(witness_label, witness_ind));
    elseif strcmp(regtype, 'l1')
        rho = 10 * bestregval;
        x = bmSteva(x0, [], [], y, ve, C, Gu, Gut, n_u, bestregval, rho, nCGD, ve_max, nIter, bmWitnessInfo(witness_label, witness_ind));
    else
        error('Invalid regularization type.');
    end
    
    % Save results
    results_folder = './reconstructions/';
    if ~exist(results_folder, 'dir')
        mkdir(results_folder);
    end
    save_path = fullfile(results_folder, sprintf('monalisa_%s_%s.mat', prefix, regtype));
    save(save_path, 'x', 'bestregval', 'prefix', 'regtype');
    fprintf('Reconstruction saved to %s\n', save_path);
end
