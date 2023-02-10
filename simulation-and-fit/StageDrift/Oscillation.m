classdef Oscillation < StageDrift

    properties
        amplitude (1,1) = 1
        frequency (1,1) = 1
        direction (1,1) {mustBeInFullRadialRange} = 0
        unit string = 'nm'
    end

    methods
        function obj = Oscillation(amplitude, frequency, direction, nSteps, unit)
            if nargin > 0
                obj.amplitude = amplitude;
                obj.frequency = frequency;
                obj.unit = unit;
                mustBeInteger(nSteps); mustBePositive(nSteps);
            end
            if nargin < 3
                nSteps = 10;
            end

            if numel(obj.amplitude) == 1
                amplitude = obj.amplitude.*ones(nSteps-1,1);
            elseif size(obj.amplitude,2) ~= nSteps
                error('Amplitude must be of length nSteps')
            else
                amplitude = obj.amplitude;
            end
            grid = (1:(nSteps-1))';

            rotationMatrix = [cos(direction), -sin(direction); sin(direction), cos(direction)];
            trajectory = [0 0; amplitude.*[sin(2*pi*grid/nSteps*frequency), zeros(nSteps-1,1)]]*rotationMatrix';
            obj.motion = Length(trajectory, obj.unit);
        end
    end

    methods (Static)
        function amplitude = createRandomAmplitudes(nSteps)
            amplitude = cumsum(randn(nSteps,2));
        end
    end
end
