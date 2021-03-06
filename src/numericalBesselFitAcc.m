function [fitresult, gof] = numericalBesselFitAcc(data, fitString,v)

% time = data.time(2:end);
% y = data.norm(1,2:end) - data.norm(1,2);
% y = y/y(end);
time = data.timeAx;
y = data.accumulation(v,:)./data.reservoir(v,:);
%y = data.accumulation(1,:)-data.refGrn;

[xData, yData] = prepareCurveData( time, y );

% Set up fittype and options.
ft = fittype( fitString, 'independent', 't', 'dependent', 'y' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
% opts.Lower = [0]; % D only
% opts.Upper = [data.D(1)];
% opts.StartPoint =[data.D(1)]; % D only
opts.Lower = [0 -Inf -Inf]; % D c1 c2
%opts.Upper = [data.D(1)];

% use these start points if there's already an estimate of D
%opts.StartPoint =[max(1,data.D(1)) 1 0]; % D c1 c2

% use these start points if there is no estimate of D
opts.StartPoint =[20 1 0]; % D c1 c2
% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'FRAP fit' );
h = plot( fitresult, xData, yData );
legend( h, 'y vs. time', 'FRAP fit', 'Location', 'NorthWest' );
% Label axes
xlabel time
ylabel('y')
grid on

end