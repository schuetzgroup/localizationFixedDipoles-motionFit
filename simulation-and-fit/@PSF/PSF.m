classdef PSF

    properties
        % PSF image
        image
        nPixels (1,1) {mustBeInteger, mustBePositive, mustBeOdd} = 17

        % Fluorophore
        dipole Dipole = Dipole(0,0)
        position (1,3) Length = Length([0 0 0], 'nm')
        nPhotons (1,1) {mustBeInteger, mustBeNonnegative} = 1e5
        shotNoise (1,1) logical = 0
        reducedExcitation (1,1) logical = 0
        stageDrift StageDrift = NoStageDrift() % relative lateral motion from start point (3rd component = defocus) 
        
        % Microscope setup
        wavelength (1,1) Length {mustBePositive} = Length(680,'nm')
        defocus (1,1) Length = Length(0, 'nm') % defocus of objective lens (negative values = moving focus into fluid)
        astigmatism = 0 % Zernike coefficient (in units of wavelength, i.e. 0.11 for Zernike coefficient 0.11*lambda)
        objectiveNA (1,1) {mustBePositive} = 0.7
        objectiveFocalLength (1,1) Length {mustBePositive} = Length(3,'mm')
        refractiveIndices (1,3) {mustBePositive} = [1.33 1.46 1] % [RI_specimen, RI_intermed, RI_immoil]
        heightIntermediateLayer (1,1) Length {mustBeNonnegative} = Length(0, 'mm')

        % Back focal plane
        backFocalPlane
        phaseMask = @(n) EmptyPhaseMask(n)
        attenuation = @(n) Aperture(n)
        nDiscretizationBFP (1,1) {mustBeInteger, mustBePositive, mustBeOdd} = 129
        
        % Camera
        pixelSize (1,1) Length {mustBePositive} = Length(108,'nm')
        pixelSensitivityMask {mustBeNonnegative, mustBeOddSize} = PixelSensitivity.uniform(9)
        backgroundNoise (1,1) {mustBeNonnegative} = 0
    end

    properties (Hidden)
        phaseMaskObj
        attenuationMaskObj
    end
    
    properties (Hidden, Access=protected)
        Defocus
        pupilMask
        fieldBFP
        chirpZTransform
    end
    properties (Dependent)
        oversampling (1,1)
        positionInPixelFromCenter (1,3)
        positionInPixelFromOrigin (1,3)
        positionInNanometerFromCenter (1,3)
    end
    properties (Dependent, Hidden)
        unitObjectSpace (1,1)
        unitKSpace (1,1)
    end

    methods
        %% Constructor
        function obj = PSF(par)
            if nargin > 0
                obj = setInputParameters('PSF', obj, par);
            end
            obj = obj.createImage();
        end

        %% Get methods for dependent properties
        function oversampling = get.oversampling(obj)
            oversampling = size(obj.pixelSensitivityMask, 1);
        end
        function positionInPixelFromCenter = get.positionInPixelFromCenter(obj)
            positionInPixelFromCenter = obj.position.inPixels(obj.pixelSize);
        end
        function positionInPixelFromOrigin = get.positionInPixelFromOrigin(obj)
            pos = obj.positionInPixelFromCenter;
            positionInPixelFromOrigin = [(obj.nPixels+1)/2 + pos(1), (obj.nPixels+1)/2 + pos(2)];
        end
        function positionInNanometerFromCenter = get.positionInNanometerFromCenter(obj)
            positionInNanometerFromCenter = obj.positionInPixelFromCenter .* obj.pixelSize.inMeter .* 1e9;
        end
        function unitObjectSpace = get.unitObjectSpace(obj)
            unitObjectSpace = obj.pixelSize.inMeter / obj.oversampling;
        end
        function unitKSpace = get.unitKSpace(obj)
            % Largest spatial frequency passed by objective lens
            maxSpatialFrequency = obj.objectiveNA / obj.wavelength.inMeter;
            maxAngularFrequency = 2*pi * maxSpatialFrequency;
            % Unit in pupil space (k-space)
            unitKSpace = 2 * maxAngularFrequency / obj.nDiscretizationBFP; % full range covers 2 times the maximum frequency
        end

        %% Functions for calculation of PSF

        function obj = createImage(obj)
            [obj.Defocus, obj.pupilMask, obj.chirpZTransform, obj.phaseMaskObj, obj.attenuationMaskObj] = obj.setup();
            
            bfp = BackFocalPlane(obj);
            obj.backFocalPlane = bfp;

            % Apply phase mask
            obj.fieldBFP.x = obj.phaseMaskObj.apply(bfp.electricField.x);
            obj.fieldBFP.y = obj.phaseMaskObj.apply(bfp.electricField.y);

            % Apply attenuation mask
            obj.fieldBFP.x = obj.attenuationMaskObj.apply(obj.fieldBFP.x);
            obj.fieldBFP.y = obj.attenuationMaskObj.apply(obj.fieldBFP.y);

            psf = zeros(obj.nPixels,obj.nPixels); 
            for k=1:size(obj.stageDrift.motion,1)
                % Apply aberrations
                aberrations = getAberrations(obj,k);
                aberratedFieldBFP = applyAberrations(obj, aberrations);
                
                % Get image from BFP field
                psf = psf + getIntensitiesCamera(obj, aberratedFieldBFP)./size(obj.stageDrift.motion,1);
            end

            psf = adjustExcitation(obj, psf);
            psf = applyShotNoise(obj, psf);
            psf = addBackgroundNoise(obj, psf);
            obj.image = psf;
        end

        function [Defocus, pupilMask, chirpZTransform, phaseMaskObj, attenuationObj] = setup(obj)
            Defocus = sphericalAberrationFromRefractiveIndexMismatch(obj);
            parMask.nGrid = obj.nDiscretizationBFP;
            pupilMask = Mask(parMask);
            chirpZTransform = ChirpZTransform(obj);
            phaseMaskObj = obj.phaseMask(obj.nDiscretizationBFP);
            attenuationObj = obj.attenuation(obj.nDiscretizationBFP);
        end

        function aberrations = getAberrations(obj,k)
            zernikeConstant = obj.unitObjectSpace * obj.unitKSpace * obj.nDiscretizationBFP /4;
            pos = obj.position.inPixels(obj.pixelSize);
            relativeMotion = obj.stageDrift.motion.inPixels(obj.pixelSize);
            
            if size(obj.stageDrift.motion,2)==3
                axialMotionInMeter = obj.stageDrift.motion.inMeter;
                axialMotionInMeter = axialMotionInMeter(:,3);
            else
                axialMotionInMeter = zeros(size(obj.stageDrift.motion,1),1); 
            end
           
            % (circle) Zernike polynomials
            Z = ZernikePolynomials.getInstance(obj.pupilMask);
            coeffPos = - (pos(1:2) + relativeMotion(k,1:2)) .* obj.oversampling .* zernikeConstant / (2*pi); % x/y-tilt
            if obj.astigmatism ~= 0
                coeffAstigmatism = obj.astigmatism; % vertical astigmatism
                aberrations = Z.getAberration([2,3,6],[coeffPos,coeffAstigmatism]) + 1/(2*pi)*(axialMotionInMeter(k)+obj.defocus.inMeter)*obj.Defocus;
            else
                aberrations = Z.getAberration([2,3],coeffPos) + 1/(2*pi)*(axialMotionInMeter(k)+obj.defocus.inMeter)*obj.Defocus;
            end
        end

        function fieldBFP = applyAberrations(obj, aberrations)
            mask = obj.pupilMask.values .* exp( 1i*2*pi*aberrations) ;%+ 1i*(obj.defocus.inMeter)*obj.Defocus );
            fieldBFP.x = obj.fieldBFP.x .* mask;
            fieldBFP.y = obj.fieldBFP.y .* mask;
        end

        function psf = adjustExcitation(obj, psf)
            % Dipole excitation probability dependent on
            % angle between electric field and dipole orientation
            if obj.reducedExcitation
                psf = psf * (sin(obj.dipole.inclination))^2;
            end
        end
        function psf = applyShotNoise(obj, psf)
            if obj.shotNoise
                psf = poissrnd(psf);
            end
        end
        function psf = addBackgroundNoise(obj, psf)
            if obj.backgroundNoise ~= 0
                psf = psf + poissrnd(obj.backgroundNoise*ones(size(psf)));
            end
        end

        % Additional functions defined in separate files
        [SA_out,Defocus,a] = sphericalAberrationFromRefractiveIndexMismatch(obj, removeDefocus)
        [psf, I_xx, I_yx] = getIntensitiesCamera(obj, mask)
        
        par = readParameters(obj)
        
        
        %% Plot function
        function plot(obj,setLimitsFullRange)
            if nargin==2 && setLimitsFullRange
                imagesc(obj.image(:,:),[0 max(obj.image(:))]);
            else
                imagesc(flip(obj.image(:,:)));
            end
            colorbar;
            axis equal; axis tight; axis xy;
            title(['Pixelsize = ', num2str(obj.pixelSize.inNanometer), 'nm. ', ...
                'Dipole = [', num2str(obj.dipole.inclination),',',num2str(obj.dipole.azimuth) ,']']);
        end
    end
end

