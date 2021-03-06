function [fitresult, gof] = FRAPfit2component(time, y)
%CREATEFIT1(TIME,BLEACH2)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : time
%      Y Output: bleach2
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 14-Jan-2019 15:04:51


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( time, y );

% Set up fittype and options.
ft = fittype( 'a1*(1-exp(-t/tau1))+a2*(1-exp(-t/tau2))+c', 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0 -Inf 0 0];
opts.StartPoint = [0.1 0.01 0.1 500 50];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'FRAP fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. time', 'FRAP fit', 'Location', 'NorthWest' );
% Label axes
xlabel time
ylabel y
grid on


