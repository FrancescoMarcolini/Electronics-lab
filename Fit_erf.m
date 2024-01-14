% This is a script for fitting the CDF (Cumulative Distribution Function):
%
%  Create a fit.
%
%      X Input : x_1Kohm...
%      Y Output: Vvolt_1Kohm...
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.
% CDF = p0/2*(1 +-  erf(sqrt(2)*(x-x0)/w0)) or 
% CDF = y0 +-  p0/2*erf(sqrt(2)*(x-x0)/w0) 
%
% x is the unknown variable
% p0,y0,w0,x0 are the fitting parameters
%
% Antonio Perreca 12/04/2017
% -------------------------------------------------------------------------
clear all
% Loading data
% Measurement made with the Thorlabs power meter (x_mm,PuW,x_mm1,PuW1)
%load('PD_S121C_potenza.mat')
% Measurement made with the Thorlabs power meter (x_mm,PuW,x_mm1,PuW1)
%load('FS100_R_load_1KOhm_500Ohm_dlaser_lama_100cm.mat')

Dati = load('datiKEM1.m');
X = Dati(:,1)-Dati(1,1);
V = Dati(:,2);

%% Fit: 'untitled fit 1'.
%[xData, yData] = prepareCurveData( x_1Kohm, Vvolt_1Kohm );
[xData, yData] = prepareCurveData( X, V );

% Set up fittype and options.
%ft = fittype( 'y0 - p0/2 * erf(sqrt(2)*(x-x0)/w0)', 'independent', 'x', 'dependent', 'z' );
ft = fittype( 'p0/2 *(1 - erf(sqrt(2)*(x-x0)/w))', 'independent', 'x', 'dependent', 'z' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Algorithm = 'Levenberg-Marquardt';
opts.Display = 'Off';

%opts.StartPoint = [0.853031117721894 0.622055131485066 0.350952380892271 0.513249539867053];
%opts.StartPoint = [10 0.350952380892271 0.513249539867053];

%opts.Display = 'Off';
%opts.StartPoint = [0.546881519204984 0.957506835434298 0.964888535199277 0.157613081677548];
%opts.StartPoint = [0.0758542895630636 0.0539501186666071 0.530797553008973 0.779167230102011];

%opts.StartPoint = [0.4 0.9 0.1 0.1];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'fit power' );
dy=ones(size(yData))*10^(-3)/sqrt(12);
h = plot( fitresult);
hold on
errorbar(xData, yData, dy,'LineStyle', 'none', 'color', 'b','CapSize',1, 'marker','.','MarkerSize', 10);
legend( 'fit', 'PuW ', 'Location', 'NorthEast' );
% Label axes
xlabel('x section [mm]')
ylabel('power [uW]')
grid on
fitresult

[wh zh]=Fit_beam_size_Trento('KemFitData.txt')
%zoom 
title('ZOOM sui dati - Gaussian Beam Profile');
xlim([0 2]);
ylim([1.6 1.9]*10^-3)
%prova per errori
outfit=Fit_ERF_function(