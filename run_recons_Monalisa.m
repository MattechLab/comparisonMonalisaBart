% generate the simulated rawdata for the comparison.
clear;
close all;
% add monalisa to path
addpath(genpath('/Users/mauroleidi/Desktop/monalisa/src')); 

regvals = [logspace(-1.5, 1.5, 40)];
prefix = 'phantom';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 0.91525
% Best regval new = 0.37751

regvals = [logspace(-1.5, 1.5, 40)];
prefix = 'phantom';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 0.76668
% Best regval new = 0.64223

regvals = [logspace(-2, 1, 40)];
prefix = 'eye';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 0.28943
% Best regval new = 0.24245

regvals = [logspace(-1, 1, 40)];
prefix = 'eye';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 1.3434
% Best regval new = 1.1938

regvals = [logspace(-3, 0, 40)];
prefix = 'cardiac';
regtype = 'l1';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 0.02309
% Best regval new = 0.028943

regvals = [logspace(-2, 0, 40)];
prefix = 'cardiac';
regtype = 'l2';
lineSearchrecon(prefix, regtype, regvals); % Best regval = 0.27283
% Best regval new = 0.43755

%% CONVERGENCE CHECK : %% MAKE CONVERGENCE PLOT FROM WITNESS_INFO

prefix = 'cardiac';
regtype = 'l2';
convergencePlot(prefix,regtype) % this automatically load the witnessinfo object

%% RUN FINAL RECONSTRUCTIONS AND SAVE RESULTS

prefix = 'phantom';
regtype = 'l1';
bestregval = 0.37751;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)

prefix = 'phantom';
regtype = 'l2';
bestregval = 0.64223;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)

prefix = 'eye';
regtype = 'l1';
bestregval = 0.24245;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)

prefix = 'eye';
regtype = 'l2';
bestregval = 1.1938;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)

prefix = 'cardiac';
regtype = 'l1';
bestregval = 0.028943;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)

prefix = 'cardiac';
regtype = 'l2';
bestregval = 0.43755;
runFinalReconstruction(prefix, regtype, bestregval)
convergencePlot(prefix,regtype)