% generate the simulated rawdata for the comparison.
clear;
close all;

regvals = [logspace(0, 2, 10)];
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


%regvals = [logspace(-4, 4, 40)];
%prefix = 'cardiac';
%regtype = 'l1';
%lineSearchrecon(prefix, regtype, regvals);


%regvals = [logspace(-4, 4, 40)];
%prefix = 'cardiac';
%regtype = 'l2';
%lineSearchrecon(prefix, regtype, regvals);

%regvals = [logspace(-4, 4, 40)];
%prefix = 'cardiac2';
%regtype = 'l1';
%lineSearchrecon(prefix, regtype, regvals);


%regvals = [logspace(-4, 4, 15)];
%prefix = 'cardiac2';
%regtype = 'l2';
%lineSearchrecon(prefix, regtype, regvals); 

regvals = [logspace(-4, 4, 40)];
prefix = 'cardiacnew';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals); %0.074438


regvals = [logspace(-4, 4, 15)];
prefix = 'cardiacnew';
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

%% RUN FINAL RECONSTRUCTIONS AND SAVE RESULTS

prefix = 'phantom';
regtype = 'l1';
bestregval = 0.651245;
runFinalReconstruction(prefix, regtype, bestregval)

prefix = 'phantom';
regtype = 'l2';
bestregval = 0.78965;
runFinalReconstruction(prefix, regtype, bestregval)

prefix = 'brain';
regtype = 'l1';
bestregval = 0.307;
runFinalReconstruction(prefix, regtype, bestregval)

prefix = 'brain';
regtype = 'l2';
bestregval = 0.046416;
runFinalReconstruction(prefix, regtype, bestregval)

%prefix = 'cardiac';
%regtype = 'l1';
%bestregval = 0.0017013;
%runFinalReconstruction(prefix, regtype, bestregval)

%prefix = 'cardiac';
%regtype = 'l2';
%bestregval = 0.30703;
%runFinalReconstruction(prefix, regtype, bestregval)


%prefix = 'cardiac2';
%regtype = 'l1';
%bestregval = 0.19145;
%runFinalReconstruction(prefix, regtype, bestregval)


%prefix = 'cardiac2';
%regtype = 'l2';
%bestregval = 0.49239;
%runFinalReconstruction(prefix, regtype, bestregval)

prefix = 'cardiacnew';
regtype = 'l1';
bestregval = 0.074438;
runFinalReconstruction(prefix, regtype, bestregval)


prefix = 'cardiacnew';
regtype = 'l2';
bestregval = 0.49239;
runFinalReconstruction(prefix, regtype, bestregval)

