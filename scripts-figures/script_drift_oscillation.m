%% Stage Drift Linear Drift 
% Calculates precision for various speeds


rng('default')

%% Specify Parameters

samplingRate = 25; 

% For saving
fullPathSave = createResultsFolder(['results/oscillation_',num2str(samplingRate)]);

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

maxAmplitude = 400;
nAmplitudePoints = 11; % (including 0)
direction = 0;  
amplitude = linspace(0,maxAmplitude,nAmplitudePoints);
frequency = 10; % amount of oscillations in one frame 

% Amount of simulations
nSimulations = 1000;


%% Run simulations

accuracy = zeros(nAmplitudePoints,3);
precision = zeros(nAmplitudePoints,3);
data = zeros(nSimulations, 3, nAmplitudePoints);

for amplitudeIndex = 1:nAmplitudePoints
    estimates = zeros(nSimulations,3);
    for k=1:nSimulations
        fprintf('Run %d...\n', k)
        parSimulation = par;

        %% Create PSF
        % Generate high resolution motion
        parSimulation.stageDrift = Oscillation(amplitude(amplitudeIndex),frequency,direction,motionSteps*resolutionFactor,'nm');

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

    precision(amplitudeIndex,:) = std(estimates);
    accuracy(amplitudeIndex,:) = mean(estimates);
    data(:,:,amplitudeIndex) = estimates; 
end


%% Save results
saveAsDat(fullPathSave, 'oscillationTest_precision', [amplitude', precision], {'amplitude','x','y','d'})
saveAsDat(fullPathSave, 'oscillationTest_accuracy', [amplitude', accuracy], {'amplitude','x','y','d'})
save(fullfile(fullPathSave, 'data_oscillationTest.mat'), 'data')
