classdef BrownianMotion < StageDrift

    properties
        diffusionConstant (1,1) = 10
        unit string = 'nm'
    end

    methods
        function obj = BrownianMotion(diffusionConstant,nSteps,unit)
            if nargin > 0
                obj.diffusionConstant = diffusionConstant;
                obj.unit = unit;
                mustBeInteger(nSteps); mustBePositive(nSteps);
            end
            if nargin < 2
                nSteps = 10;
            end
            
            % Brownian motion corresponds to normally distributed displacements
            % Stepwidth:
            % sigma^2 = 2*dim*D*t
            % = 4*D*t (for 2-dimensional diffusion)
            sigma = sqrt(4*obj.diffusionConstant);
            steps = normrnd(0, sigma, nSteps-1, 2);
            trajectory = [0 0; cumsum(steps)];
            obj.motion = Length(trajectory, obj.unit);
        end
    end
end