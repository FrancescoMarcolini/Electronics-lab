function setbarwidth(h, x, xlim)
  % SETBARWIDTH  Set width of errorbars in plots produced by ERRORBAR.
  %
  % SETBARWIDTH(HANDLE, XVALS, WIDTH)
  %
  %   HANDLE - errorbar plot handle as returned by ERRORBAR
  %   XVALS  - x values of the data points in the plot
  %   WIDTH  - half width of errobars. can be specified as a single
  %            scalar or as a vector with the same lenght as XVALS.
  %
  % See also ERRORBAR.
  %
  % Code courtesy of D. Nicolodi and Matlab Central
  %
  % $Id$
  
  if verLessThan('matlab', '8.4')
    % Matlab R2012b
    h0 = get(h, 'Children');
    b = get(h0(2), 'Xdata');
    
    switch length(xlim)
      case 1
        xlim = xlim * ones(size(x));
      case length(x)
        % nothing to do
      otherwise
        error('WIDTH should be a scalar or have the same size as XVALS')
    end
    
    for kk = 1:length(x)
      temp_line = b(9*(kk-1)+1 : 9*(kk));
      temp_line(4) = x(kk) - xlim(kk);
      temp_line(5) = x(kk) + xlim(kk);
      temp_line(7) = x(kk) - xlim(kk);
      temp_line(8) = x(kk) + xlim(kk);
      b(9*(kk-1)+1 : 9*(kk)) = temp_line;
    end
    
    assert(all(size(get(h0(2), 'Ydata')) == size(b)));
    
    set(h0(2), 'Xdata', b)
    
  else
    % Matlab R2014b
    switch length(xlim)
      case 1
        xlim = xlim * ones(2 * numel(x), 1);
      case length(x)
        % nothing to do
        for jj = 1:numel(xlim)
          xlim_n = [xlim_n xlim(jj) xlim(jj)];
        end
        xlim = xlim_n;
      otherwise
        error('WIDTH should be a scalar or have the same size as XVALS')
    end
    
    % Check that the drawing actually happens
    kk = 0;
    while isempty(h.Bar.VertexData) && kk < 100
      pause(0.05);
      kk = kk + 1;
    end
    b = h.Bar.VertexData(1, :);
    
    % the error bar 'caps' coordinats come after a list of 2 data points per error bar
    N = 2*length(x);
    
    m = 1; % Multiplier for addition/subtraction
    for ii = 1:N
      m = -1*m; % Switch between subtraction and addition
      b(N+ii:N:end) = b(ii) + m*xlim(ii); % Change xdata with respect to the chosen ratio
    end
    h.Bar.VertexData(1, :) = b;
  end
