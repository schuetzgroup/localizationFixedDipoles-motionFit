%% Script for creating plots

close all
clear all

%% Static PSFs

par.dipole = Dipole(pi/2,0);

figure
subplot(2,4,1) 
psf = PSF(par); 
limit = [0 max(psf.image(:))];
imagesc(psf.image, limit)
title('\theta = \pi/2', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

par.dipole = Dipole(pi/3,0); 
subplot(2,4,2) 
psf = PSF(par); 
imagesc(psf.image, limit)
title('\theta = \pi/3', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square

par.dipole = Dipole(pi/6,0); 
subplot(2,4,3) 
psf = PSF(par); 
imagesc(psf.image, limit)
title('\theta = \pi/6', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

par.dipole = Dipole(0,0); 
subplot(2,4,4) 
psf = PSF(par); 
imagesc(psf.image, limit)
title('\theta = 0', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

 
% With stage drift

clear par

par.stageDrift = LinearDrift(400/500,pi/4, 500,'nm');
X = double(par.stageDrift.motion); 

par.dipole = Dipole(pi/2,0); 
subplot(2,4,5) 
psf = PSF(par); 
imagesc(psf.image, limit)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

par.dipole = Dipole(pi/3,0); 
subplot(2,4,6) 
psf = PSF(par); 
imagesc(psf.image, limit)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

par.dipole = Dipole(pi/6,0); 
subplot(2,4,7) 
psf = PSF(par); 
imagesc(psf.image, limit)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

par.dipole = Dipole(0,0); 
subplot(2,4,8) 
psf = PSF(par); 
imagesc(psf.image, limit)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square


%% Stage drift

clear par 

par.dipole = Dipole(pi/2,0);

figure
subplot(1,4,1) 
psf = PSF(par); 
imagesc(psf.image)
title('Static', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

subplot(1,4,2) 
par.stageDrift = LinearDrift(400/500,pi/4, 500,'nm');
psf = PSF(par); 
X = double(par.stageDrift.motion); 
imagesc(psf.image)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
title('Linear Drift', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

subplot(1,4,3) 
rng(13)
par.stageDrift = BrownianMotion(45, 500,'nm');
psf = PSF(par); 
X = double(par.stageDrift.motion); 
imagesc(psf.image)
hold on
plot(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'Color','black')
title('Diffusion', 'FontSize', 15)
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square 

subplot(1,4,4) 
rng(14)
par.stageDrift = Oscillation(300,10,0, 500,'nm');
psf = PSF(par); 
X = double(par.stageDrift.motion); 
imagesc(psf.image)
title('Oscillation', 'FontSize', 15)
hold on
scatter(9+10^7.*X(:,1), 9-10^7.*X(:,2), 'black','.')
set(gca,'XTickLabel',[]);
set(gca,'YTickLabel',[]);
axis square
