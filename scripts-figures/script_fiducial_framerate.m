%% Stage Drift Brownian Motion
% Calculates precision for various diffusion constants

rng('default')

%% Specify Parameters
RMSD = 100; 

% For saving
fullPathSave = createResultsFolder(['results/fiducialFrameRate_', num2str(RMSD)]);

% Fluorophore
par.shotNoise = 1;

% Microscope setup
par.astigmatism = 0.11;

% Camera
par.nPixels = 17;
par.backgroundNoise = 100;

% Stage drift (Brownian motion) 
resolutionFactor = 20;
FiducialSamplingPoints = [1;2;4;5;10;20;25]; 

nFiducialSamplingPoints = length(FiducialSamplingPoints);
% Amount of simulations
nSimulations = 2000;


%% Run simulations
% Create arrays for storing results 
% localization accuracy
accuracy = zeros(nFiducialSamplingPoints,3);
% localization precision
precision = zeros(nFiducialSamplingPoints,3);
% errors for each individual simulation
data = zeros(nSimulations,3,nFiducialSamplingPoints);

for fiducialIndex = 1:nFiducialSamplingPoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)

        motionSteps = FiducialSamplingPoints(fiducialIndex); 
        samplingRate = motionSteps; 

        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        parSimulation.stageDrift = BrownianMotion((RMSD/2).^2/(motionSteps*resolutionFactor),motionSteps*resolutionFactor,'nm');

        % Generate random dipole orientation, position and defocus
        parFluorophore = createRandomFluorophore();
        parSimulation = appendStruct(parSimulation, parFluorophore);

        psf = PSF(parSimulation);

        %% Fit PSF
        % Pass drift
        parFit.stageDrift = parSimulation.stageDrift.reduceSamplingRate(samplingRate);
        parFit.stageDrift = parFit.stageDrift.addNoise(sqrt(samplingRate)); 

        % Create angle and background noise estimates
        withNoise = 1;
        parFitEstimates = createParameterEstimates(parSimulation, withNoise);
        parFit = appendStruct(parFit, parFitEstimates);

        est = FitPSF(psf,parFit);

        %% Evaluate error
        pos = parSimulation.position.inNanometer;
        estimates(k,:) = [est.estimatesPositionDefocus.ML] - [pos(1) pos(2) parSimulation.defocus.inNanometer];
    end

    precision(fiducialIndex,:) = std(estimates);
    accuracy(fiducialIndex,:) = mean(estimates);
    data(:,:,fiducialIndex) = estimates; 
end


%% Save results

saveAsDat(fullPathSave, 'fiducialFrameRateTest_precision', [FiducialSamplingPoints, precision], {'Sampling','x','y','d'})
saveAsDat(fullPathSave, 'fiducialFrameRateTest_accuracy', [FiducialSamplingPoints, accuracy], {'Sampling','x','y','d'})
save(fullfile(fullPathSave,'data_fiducialFrameRateTest.mat'), 'data')

