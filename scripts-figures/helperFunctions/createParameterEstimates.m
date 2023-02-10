function estimates = createParameterEstimates(parPsf, withNoise)
    arguments
        parPsf struct
        withNoise logical = 0
    end
    if withNoise
        estimates.noiseEstimate = max([mean(poissrnd(parPsf.backgroundNoise*ones(parPsf.nPixels)),'all'), 1e-5]);
        estimates.angleInclinationEstimate = addAngleEstimateError(parPsf.dipole.inclination, 2);
        estimates.angleAzimuthEstimate = addAngleEstimateError(parPsf.dipole.azimuth, 2);
    else
        estimates.noiseEstimate = max([parPsf.backgroundNoise,1e-5]);
        estimates.angleInclinationEstimate = parPsf.dipole.inclination;
        estimates.angleAzimuthEstimate = parPsf.dipole.azimuth;
    end
end