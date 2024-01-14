


function [wh zh]=Fit_beam_size_Trento(filename)
%function [wh zh wv zv]=Fit_beam_size_Trento(filename)

%
% This is a fit function of a gaussian beam provided by two sets of data:
% Horizontal and vertical sections.
%
% Required Functions:
% * [xv0_new xvr_new]=opt_st_point(x,wv) 
%
% INPUT:
% * 'filename.txt/m' which is a file.txt where the data are stored;
% 
% OUTPUT:
% * wh = Horizontal beam size
% * zh = Horizontal waist location
% * Wv = Vertical beam size
% * zv = Vertical waist location
%
%
% Antonio Perreca 23-02-2013
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%   load data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

s=load (filename,'-ascii')

%the factor 0.0254 has to be used only if the distances are in inches
d=s(:,1)*0.01;%*0.0254;
% The values are converted in meters 
%wh=(s(:,2)/2)*1e-6;
%wv=(s(:,3)/2)*1e-6;
wh=(s(:,2))*1e-3;% when the measurement is in mm
%wv=(s(:,3))*1e-3;

%%%%%%%%%%%%%%%%%   Fitting data   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%axes1 = axes('Parent',figure,'FontSize',24);
% function [w0,z0] = FT_fit_beam_size2(lambda, w0, z0, w_data, z_data, tolerance, show_plot,color_data,color_fit)
[wh,zh] = FT_fit_beam_size2(532e-9, 800e-6, 1, wh, d, 1e-6, 1,'ro','b-');
legend('x-section data',['x-fit: w0 = ' num2str(wh),' z0 = ' num2str(zh) ])

% hold on
% function [w0,z0] = FT_fit_beam_size2(lambda, w0, z0, w_data, z_data, tolerance, show_plot,color_data,color_fit)
% [wv,zv] = FT_fit_beam_size2(520e-9, 800e-6, 1, wv, d, 1e-6, 1,'mo','r-')
% legend('x-section data',['x-fit: w0 = ' num2str(wh),' z0 = ' num2str(zh) ],'y-section data',['y-fit: w0 = ' num2str(wv),' z0 = ' num2str(zv) ])

title('Gaussian beam profile')
%%%%%%%%%%%%%%%%%%   Print results   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
%clc
  disp('--------------------------- Beam parameters -----------------------------------')
  
  disp(' ')
  
    disp(['Horizontal waist size = ' num2str(wh) ' m']);
    disp(['Horizontal beam position = ' num2str(zh) ' m']);
    %disp(['Vertical waist size = ' num2str(wv) ' m']);
    %disp(['Vertical beam position = ' num2str(zv) ' m']);
   % disp(['Lamp_Voltage_' num2str(jj) ' = ' num2str(meanUV(jj).data.y) ' +/- '  num2str(meanUV(jj).data.dy) ' V'])
  
  disp(' ')
  
  disp('-------------------------------------------------------------------------------')
  