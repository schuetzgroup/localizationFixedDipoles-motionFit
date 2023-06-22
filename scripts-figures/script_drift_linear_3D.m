%% Stage Drift Linear Drift 
% Calculates precision for various speeds


rng('default')

%% Specify Parameters
% sampling rate
samplingRate = 25; 
% For saving
fullPathSave = createResultsFolder(['results/linearDrift3D_',num2str(samplingRate)]);

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

direction = 0;  
speed = 200;


nAxialDriftPoints = 11; % (including 0)
axialSpeed = linspace(0,400,nAxialDriftPoints);

% Amount of simulations
nSimulations = 1000;


%% Run simulations
% Create arrays for storing results 
% localization accuracy
accuracy = zeros(nAxialDriftPoints,3);
% localization precision
precision = zeros(nAxialDriftPoints,3);
% errors for each individual simulation
data = zeros(nSimulations,3,nAxialDriftPoints);

for speedIndex = 1:nAxialDriftPoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)
        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        lateralDrift = LinearDrift(speed/(motionSteps*resolutionFactor),direction,motionSteps*resolutionFactor,'nm');

        % generate axial oscillation (in nanometer) 
        axialDrift = Length(axialSpeed(speedIndex)*(linspace(0, 1, motionSteps*resolutionFactor))', 'nm'); 


        parSimulation.stageDrift = StageDrift(Length([lateralDrift.motion.inNanometer, axialDrift.inNanometer],'nm')); 

        % Generate random dipole orientation, position and defocus
        parFluorophore = createRandomFluorophore();
        parSimulation = appendStruct(parSimulation, parFluorophore);

        psf = PSF(parSimulation);

        %% Fit PSF
        % Pass drift
        parFit.stageDrift = parSimulation.stageDrift.deleteAxialDrift();
        parFit.stageDrift = parFit.stageDrift.reduceSamplingRate(samplingRate);
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

    precision(speedIndex,:) = std(estimates);
    accuracy(speedIndex,:) = mean(estimates);
    data(:,:,speedIndex) = estimates; 
end


%% Save results
saveAsDat(fullPathSave, 'linearDriftTest_precision', [amplitude', precision], {'speed','x','y','d'})
saveAsDat(fullPathSave, 'linearDriftTest_accuracy', [amplitude', accuracy], {'speed','x','y','d'})
save(fullfile(fullPathSave,'data_linearDrift.mat'), 'data')
