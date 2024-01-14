%
%---------------------------------------------------------------------  
% function [w0,z0] = FT_fit_beam_size2(lambda, w0, z0, w_data, z_data, tolerance, show_plot,color_data,color_fit)
% 
% Matlab function that fits the beam size function w(z) to a
% set of measured beam sizes of a Gaussian beam to determine
% the beam waist size and position.
%
% lambda (real):       wavelength of laser light [m]
% w0 (real):           waist radius (x-z,y-z plane) [m]
% z0 (real):           position on z-axis (distance from beam waist at z=0) [m]                    
%                      negative values: running towards the waist
%                      positive values: running away from the waist
%
%   Both w0 and z0 must be given as an initial guess to the fit and the
%   function returns the updates values after the fit.
%
% w_data (real):       vector of measured beam radi [m]
% z_data (real):       vector of positions at which w_data was measured [m]
% tolerance (real):    tolerance for Matlab function fminsearch (use e.g. 1e-6)
% show_plot (integer): if showplot=1 the function plots the data and the fittet
%                      result
% color_data ('ro'...) color of data
% color_fit ('b-'...)  color of the fitting curve
%
% Part of the SimTools package
% Andreas Freise, 05.01.2010 afreise@googlemail.com
% Updated by Antonio Perreca, 23.02.2013 antonio.perreca@hotmail.co.uk
%---------------------------------------------------------------------  
%

% Description: Fits the beam size function to a set of measured sizes of a Gaussian
% Description: beam to determine the beam waist size and position.
% Keywords: size, Gaussian, waist, position


function [w0,z0] = FT_fit_beam_size2(lambda,w0,z0,w_data, z_data, tolerance, show_plot,color_data,color_fit)
  
  
  params(1)=w0;
  params(2)=z0;

  options=optimset('Display','off', 'TolX', tolerance, 'TolFun',tolerance, 'MaxIter', 1000);
  
  [paramsout,fval,exitflag,output]=fminsearch(@int_w_z,params,options,lambda, w_data, z_data);
    
  residual=int_w_z(paramsout,lambda, w_data,z_data)
  
  fprintf(' -- Fit completed after %d iterations--\n',output.iterations);
  fprintf(' Started with: w0=%g, z=%g\n',w0,z0);

  w0=paramsout(1)
  z0=paramsout(2)

  fprintf(' Fit results:  w0=%g, z=%g\n',w0,z0);
  fprintf(' Residual = %g\n',residual);

  if (show_plot==1)
    %figure1=figure()
    %axes1 = axes('Parent',figure,'FontSize',24);
    x1=z_data;
    y1=w_data;
    
    xlimits = 15.5;
    %x2=linspace(min(z_data),max(z_data),200);
    x2=linspace(z0-xlimits,z0+xlimits,300);
    zr=pi*w0^2/lambda;
    y2=w0.*sqrt(1+((x2-z0)/zr).^2);
    hold on;
    plot(x1,y1,color_data,'LineWidth',2); %'ro'
    plot(x2,y2,color_fit,'LineWidth',2);
    grid on
    ylabel('Beam radius [m]');
    xlabel('Distance [m]');
    xlim([z0-xlimits z0+xlimits]);
    hold on
  end
  
hold on
  
  function delta_w = int_w_z(params, lambda, w_data, z_data)
    w0=params(1);
    z0=params(2);
    zr=pi*w0^2/lambda;
    w=w0.*sqrt(1+((z_data-z0)/zr).^2);
    delta_w=sqrt(sum((w_data-w).^2));
    
    
    
    
