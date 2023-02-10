classdef LinearDrift < StageDrift

    properties
        speed (1,1) = 1
        direction (1,1) = 0
        unit string = 'nm'
    end

    methods
        function obj = LinearDrift(speed, direction, nSteps, unit)
            if nargin > 0
                obj.speed = speed;
                obj.direction = direction;
                obj.unit = unit;
                mustBeInteger(nSteps); mustBePositive(nSteps);
            end
            if nargin < 3
                nSteps = 10;
            end
            
            grid = 1:(nSteps-1);
            rotationMatrix = [cos(direction), -sin(direction); sin(direction), cos(direction)];
            trajectory = [0 0; obj.speed.*grid', zeros(nSteps-1,1)]*rotationMatrix';
            obj.motion = Length(trajectory, obj.unit);
        end
    end
end