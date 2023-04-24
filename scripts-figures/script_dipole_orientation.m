%% Stage Drift Linear Drift 
% Calculates precision for various speeds


rng('default')

%% Specify Parameters
% Fluorophore
inclination = pi/3; % dipole inclination angle 
par.shotNoise = 1;

% For saving
fullPathSave = createResultsFolder(['results/orientationTest_',num2str(inclination)]);


% Microscope setup
par.astigmatism = 0.11;

% Camera
par.nPixels = 17;
par.backgroundNoise = 100;

% Stage drift
resolutionFactor = 20;
motionSteps = 25; % amount of motion steps
samplingRate = motionSteps; 

maxSpeed = 600;
nSpeedPoints = 11; % (including 0)
direction = 0;  
speed = linspace(0,maxSpeed,nSpeedPoints);

% Amount of simulations
nSimulations = 1000;


%% Run simulations
% Create arrays for storing results 
% localization accuracy
accuracy = zeros(nSpeedPoints,3);
% localization precision
precision = zeros(nSpeedPoints,3);
% errors for each individual simulation
data = zeros(nSimulations,3,nSpeedPoints);

for speedIndex = 1:nSpeedPoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)
        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        parSimulation.stageDrift = LinearDrift(speed(speedIndex)/(motionSteps*resolutionFactor),direction,motionSteps*resolutionFactor,'nm');

        % Generate random dipole orientation, position and defocus
        parFluorophore = createRandomFluorophore();
        parFluorophore.dipole.inclination = inclination; 
        parSimulation = appendStruct(parSimulation, parFluorophore);

        psf = PSF(parSimulation);

        %% Fit PSF
        % Pass static drift
        parFit.stageDrift = parSimulation.stageDrift.reduceSamplingRate(samplingRate);

        % Create angle and background noise estimates
        withNoise = 1;
        parFitEstimates = createParameterEstimates(parSimulation, withNoise);
        parFit = appendStruct(parFit, parFitEstimates);

        est = FitPSF(psf,parFit);

        %% Evaluate error
        pos = parSimulation.position.inNanometer;
        estimates(k,:) = [est.estimatesPositionDefocus.ML] - [pos(1) pos(2) parSimulation.defocus.inNanometer];
    end

    precision(speedIndex,:) = std(estimates);
    accuracy(speedIndex,:) = mean(estimates);
    data(:,:,speedIndex) = estimates; 
end


%% Save results
saveAsDat(fullPathSave, 'linearDriftTest_precision', [speed', precision], {'speed','x','y','d'})
saveAsDat(fullPathSave, 'linearDriftTest_accuracy', [speed', accuracy], {'speed','x','y','d'})
save(fullfile(fullPathSave,'data_linearDrift.mat'), 'data')
