%% finding directories
% long file name of running script
script_file = matlab.desktop.editor.getActiveFilename;

% identify directories
script_dir  = fileparts(script_file);
prev_dir    = fileparts(script_file);
prev_dir    = fileparts(prev_dir);
prev_dir    = fileparts(prev_dir);
working_dir = fileparts(prev_dir);

clear prev_dir
%% loading data
load([working_dir, filesep, 'images', filesep, 'C_3.mat']);
load([working_dir, filesep, 'images', filesep, 'h_3.mat']);

bmImage(cat(2, real(C), imag(C)))
bmImage(h)


%% accquisition and recon parameters

N       = 512; % Size of the image (e.g., 256x256)
FoV     = [300, 300];
dK_u    = 1./FoV;
N_u     = [256, 256];
n_u     = N_u;
nLines  = 30;


%% y, t, ve

t = bmTraj_fullRadial2_lineAssym2(N,nLines, dK_u(1));

nPt = size(t, 2);

y = bmSimulateMriData(h, C, t, N_u, n_u, dK_u);


ve      = bmVolumeElement(t, 'voronoi_full_radial2');
ve_max  = 3*prod(dK_u(:));
ve      = min(ve, ve_max);


%% gridded recon
x0 = bmMathilda(y, t, ve, C, N_u, n_u, dK_u);
bmImage(x0)

%% preparing gridding matrices
[Gu, Gut] = bmTraj2SparseMat(t, ve, N_u, dK_u);

%% sleva

cd([working_dir, filesep, 'recon_bastien\recon_cardiac\data\sleva_result'])
nIter = 35;
nCGD = 4;
delta = 0.5; 
witness_ind     = 1:5:nIter; % indices to save the witness
witness_label   = 'sleva_d0p5'; % label to save the file
regul_mode      = 'normal'; 

x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, delta, regul_mode, nCGD, ve_max, ...
    nIter, bmWitnessInfo(witness_label, witness_ind));


%% steva

cd([working_dir, filesep, 'recon_bastien\recon_cardiac\data\steva_result'])
nIter = 35;
nCGD = 4;
delta = 0.005; 
witness_ind     = 1:5:nIter; % indices to save the witness
witness_label   = 'sleva_d0p005'; % label to save the file
regul_mode      = 'normal'; 

rho = 10*delta; % this should be some what adapted too.
x = bmSteva(x0, [], [], y, ve, C, Gu, Gut, n_u, delta, rho, nCGD, ve_max, ...
    nIter, bmWitnessInfo(witness_label, witness_ind));
