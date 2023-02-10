function angleEstimate = addAngleEstimateError(angle, std)
    noise = normrnd(0,std)*pi/180;
    if abs(noise) > 4*std
        noise = sign(noise) * 4*std*pi/180;
    end
    angleEstimate = mod(angle + noise, 2*pi);
end