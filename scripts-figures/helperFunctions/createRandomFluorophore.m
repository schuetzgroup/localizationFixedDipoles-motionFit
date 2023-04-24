function par = createRandomFluorophore(pixelSize)
    if nargin==0
        pixelSize = 108;
    end
    par.defocus = Length(-500+1000*rand(), 'nm');
    par.dipole = Dipole(asin(rand()), 2*pi*rand());
    xpos = 2*pixelSize*rand()-pixelSize;
    ypos = 2*pixelSize*rand()-pixelSize;
    par.position = Length([xpos ypos 0], 'nm');
end