%% Stage Drift Brownian Motion
% Calculates precision for various diffusion constants
rng('default')

%% Specify Parameters

samplingRate = 25; 
% For saving
fullPathSave = createResultsFolder(['results/brownianMotion_',num2str(samplingRate)]);

% Fluorophore
par.shotNoise = 1;

% Microscope setup
par.astigmatism = 0.11;

% Camera
par.nPixels = 17;
par.backgroundNoise = 100;

% Stage drift
resolutionFactor = 20;
motionSteps = 25; % amount of motion steps
% samplingRate = motionSteps; 

maxRMSD = 400;
nDiffusionPoints = 11; % (including 0)
diffusionConstant = (linspace(0,maxRMSD/4,nDiffusionPoints)*2).^2;

% Amount of simulations
nSimulations = 1000;

%% Run simulations
% Create arrays for storing results 
% localization accuracy
accuracy = zeros(nDiffusionPoints,3);
% localization precision
precision = zeros(nDiffusionPoints,3);
% errors for each individual simulation
data = zeros(nSimulations,3,nDiffusionPoints);

for diffusionIndex = 1:nDiffusionPoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)

        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        parSimulation.stageDrift = BrownianMotion(diffusionConstant(diffusionIndex)/(motionSteps*resolutionFactor),motionSteps*resolutionFactor,'nm');

        % Generate random dipole orientation, position and defocus
        parFluorophore = createRandomFluorophore();
        parSimulation = appendStruct(parSimulation, parFluorophore);

        psf = PSF(parSimulation);

        %% Fit PSF
        % Pass static drift
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

    precision(diffusionIndex,:) = std(estimates);
    accuracy(diffusionIndex,:) = mean(estimates);
    data(:,:,diffusionIndex) = estimates; 
end


%% Save results
RMSD = sqrt(4*diffusionConstant)';
saveAsDat(fullPathSave, 'brownianMotionTest_precision', [RMSD, precision], {'RMSD','x','y','d'})
saveAsDat(fullPathSave, 'brownianMotionTest_accuracy', [RMSD, accuracy], {'RMSD','x','y','d'})
save(fullfile(fullPathSave,'data_brownianMotionTest.mat'), 'data')
