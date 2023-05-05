%% off-switching during frame is considered 
% Stage Drift Brownian Motion

rng('default')

%% Specify Parameters

% For saving
fullPathSave = createResultsFolder('results/blinking');

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

speed = 200; 
direction = 0; 

nPoints = 5;
Time = linspace(20,100,nPoints);

% Amount of simulations
nSimulations = 1000;


%% Run simulations

accuracy = zeros(nPoints,3);
precision = zeros(nPoints,3);
data = zeros(nSimulations, 3, nPoints);

for blinkIndex = 1:nPoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)

        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        parSimulation.stageDrift = LinearDrift(speed/(motionSteps*resolutionFactor),direction,motionSteps*resolutionFactor,'nm');
        parFit.stageDrift = parSimulation.stageDrift.reduceSamplingRate(motionSteps); 

        % cut off simulation stage drift after some point 
        parSimulation.stageDrift = parSimulation.stageDrift.cutOff(resolutionFactor*Time(blinkIndex)*motionSteps/100); 

        % Generate random dipole orientation, position and defocus
        parFluorophore = createRandomFluorophore();
        parSimulation = appendStruct(parSimulation, parFluorophore);

        psf = PSF(parSimulation);

        %% Fit PSF
        % Create angle and background noise estimates
        withNoise = 1;
        parFitEstimates = createParameterEstimates(parSimulation, withNoise);
        parFit = appendStruct(parFit, parFitEstimates);

        est = FitPSF(psf,parFit);

        %% Evaluate error
        pos = parSimulation.position.inNanometer;
        estimates(k,:) = [est.estimatesPositionDefocus.ML] - [pos(1) pos(2) parSimulation.defocus.inNanometer];
    end

    precision(blinkIndex,:) = std(estimates);
    accuracy(blinkIndex,:) = mean(estimates);
    data(:,:,blinkIndex) = estimates; 
end


%% Save results
saveAsDat(fullPathSave, 'blinkingTest_precision', [Time', precision], {'Time','x','y','d'})
saveAsDat(fullPathSave, 'blinkingTest_accuracy', [Time', accuracy], {'Time','x','y','d'})
save(fullfile(fullPathSave,'data_blinkingTest.mat'), 'data')

