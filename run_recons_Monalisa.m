% generate the simulated rawdata for the comparison.
clear;
close all;

regvals = [logspace(-4, 4, 40)];
prefix = 'phantom';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals);


regvals = [logspace(-4, 4, 40)];
prefix = 'phantom';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals);


regvals = [logspace(-4, 4, 40)];
prefix = 'brain';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals);

regvals = [logspace(-4, 4, 40)];
prefix = 'brain';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals);


regvals = [logspace(-4, 4, 40)];
prefix = 'cardiac';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals);


regvals = [logspace(-4, 4, 40)];
prefix = 'cardiac';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals);

regvals = [logspace(-4, 4, 40)];
prefix = 'cardiac2';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals);


regvals = [logspace(-4, 4, 40)];
prefix = 'cardiac2';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals);

%% CONVERGENCE CHECK (just an example on how you can check convergence of your reconstructions, NOT NEEDED HERE)
nIter = 35;
nCGD = 4;
witness_ind = 1:5:nIter;
frameSize = n_u;
rho = 10*delta; % this should be somewhat adapted too.
witness_info = bmWitnessInfo(witness_label, witness_ind);
witness_info.save_witnessIm_flag = true;
x = bmSleva(x0, y, ve, C, Gu, Gut, n_u, bestregval, regul_mode, nCGD, ve_max, ...
                    nIter,witness_info );
% Let's have a look
whatever = load('wit_label.mat');
bmImage(whatever.witnessInfo.witness_im)
bmImage(image - whatever.witnessInfo.witness_im)