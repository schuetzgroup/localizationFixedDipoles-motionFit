%% Example script for simulating and fitting PSF

close all
clear all

% Set input parameters
angleInclination = pi/16;
par.position = Length([0 0 0].*100,'nm');
par.dipole = Dipole(angleInclination, 0);
par.stageDrift = LinearDrift(10,1,10,'nm');
par.defocus = Length(-500+1000*rand(), 'nm');

% Simulate psf
psf = PSF(par);

% Fitting
parEst.angleInclinationEstimate = angleInclination;
parEst.angleAzimuthEstimate = 0;
parEst.stageDrift = par.stageDrift;
fitResult = FitPSF(psf, parEst);

% Plot
figure
plot(fitResult)
