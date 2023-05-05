%% MotionFit Test Script
% Calculates precision for varying levels of noise in the stageDrift estimate

rng('default')

%% Specify Parameters

samplingRate = 25; 

% For saving
fullPathSave = createResultsFolder(['results/noisyFiducialTest_',num2str(samplingRate)]);


% Fluorophore
par.shotNoise = 1; 

% Microscope setup
par.astigmatism = 0.11;

% Camera
par.nPixels = 17;
par.backgroundNoise = 100; 

% Stage drift 
speed = 400; 
direction = 0;  

resolutionFactor = 20;
motionSteps = 25; % amount of motion steps 

maxNoiseLevel = 5; % maximum total localization precision 
noiseStep = 0.5; 

nSimulations = 2000; % amount of simulations


% localization accuracy
accuracy = zeros(maxNoiseLevel/noiseStep+1,3);
% localization precision
precision = zeros(maxNoiseLevel/noiseStep+1,3);
% errors for each individual simulation
data = zeros(nSimulations,3,maxNoiseLevel/noiseStep+1);
index = 0;

for noiselevel = 0:noiseStep:maxNoiseLevel 
    estimates = zeros(nSimulations,3);
    index = index + 1; 

    for k=1:nSimulations
        fprintf('Run %d...\n', k)

        %% Simulate PSF

        parSimulation = par; 

        % Generate high resolution motion 
        parSimulation.stageDrift = LinearDrift(speed/(motionSteps*resolutionFactor),direction,motionSteps*resolutionFactor,'nm');
      
        % Generate random defocus, dipole and position 
        parFluorophore = createRandomFluorophore();
        parSimulation = appendStruct(parSimulation, parFluorophore);

        % Calculate PSF 
        psf = PSF(parSimulation); 
        

        %% Fit PSF
        parFit.stageDrift = parSimulation.stageDrift.reduceSamplingRate(samplingRate);
        parFit.stageDrift = parFit.stageDrift.addNoise(noiselevel*sqrt(samplingRate));

        % Create angle and background noise estimates
        withNoise = 1;
        parFitEstimates = createParameterEstimates(parSimulation, withNoise);
        parFit = appendStruct(parFit, parFitEstimates);

        est = FitPSF(psf,parFit);

        %% Evaluate error
        pos = parSimulation.position.inNanometer;
        estimates(k,:) = [est.estimatesPositionDefocus.ML] - [pos(1) pos(2) parSimulation.defocus.inNanometer];
    end

    precision(index,:) = std(estimates);
    accuracy(index,:) = mean(estimates);
    data(:,:,index) = estimates; 
end

%% Save results

noiselevels = (0:maxNoiseLevel/noiseStep)';
saveAsDat(fullPathSave, 'noisyFiducialTest_precision', [noiselevels, precision], {'noiselevel','x','y','d'})
saveAsDat(fullPathSave, 'noisyFiducialTest_accuracy', [noiselevels, accuracy], {'noiselevel','x','y','d'})
save(fullfile(fullPathSave,'data_noisyFiducial.mat'), 'data')

