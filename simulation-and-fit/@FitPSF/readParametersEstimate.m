function par = readParametersEstimate(psf)

    par.nPixels = psf.nPixels;
    
    % Fluorophore
    par.shotNoise =  0; 
    par.reducedExcitation = 0;

    % Microscope setup
    par.wavelength = psf.wavelength;
    par.astigmatism = psf.astigmatism;
    par.objectiveNA = psf.objectiveNA;
    par.objectiveFocalLength = psf.objectiveFocalLength;
    par.refractiveIndices = psf.refractiveIndices;
    par.heightIntermediateLayer = psf.heightIntermediateLayer;

    % Back focal plane
    par.phaseMask = psf.phaseMask;
    par.nDiscretizationBFP = psf.nDiscretizationBFP;

    % Camera
    par.pixelSize = psf.pixelSize;
end