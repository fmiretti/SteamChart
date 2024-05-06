clear
close all

% Help files:
%   help SteamChart
%   help SteamChart/add_marker
%   help SteamChart/remove_marker
%   help SteamChart/save_chart

%% Examples
% Create chart 
ch = SteamChart('sLim', [7,8], 'hLim', [2500, 3500]);

ch.isobarics = [0.06 0.14 0.4 0.8 1.4 2.6 6 12 24 46 90];
ch.isothermals = [40 80 160 200 240 280 320 380 420 460];
ch.isochorics = [];
ch.isoquality = [0.95 0.98];

% Add markers on the chart 
ch.add_marker('A', 'pT', 10, 300)
ch.add_marker('B', 'pT', 10, 200)
ch.add_marker('C', 'hs', 2725, 7.5)
ch.add_marker('D', 'ps', 1, 7.5)

% Add several markers in one function call
ch.add_marker({'D', 'E'}, 'pT', [10, 10], [220, 230])

% Remove markers
ch.remove_marker('B')
ch.remove_marker({'D', 'A'})

% Export chart to pdf
ch.save_chart('chart.pdf')

