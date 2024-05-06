classdef SteamChart < matlab.graphics.chartcontainer.ChartContainer
    % SteamChart
    %   Create printable A4 steam charts
    %   ch = SteamChart() creates a SteamChart plot and returns a handle to the
    %   chart
    %   ch = SteamChart('Name', value, ...) creates a SteamChart plot while
    %   setting the specified Name-Value properties.
    %
    %   Add markers to the chart by using the add_marker and remove_marker
    %   methods
    %
    %   Name-Value Properties
    %   -------
    %   'sLim':     1x2 numeric
    %       The entrophy axis limits
    %   'hLim':     1x2 numeric
    %       The enthalpy axis limits
    %   'sMajorGridStep', 'hMajorGridStep':     scalar
    %       The major grids spacing
    %   'sMinorGridStep', 'hMinorGridStep':     scalar
    %       The minor grids spacing
    %   'isobarics':    numeric vector
    %       Pressure values (bar) for the isobarics to draw
    %   'isothermals':    numeric vector
    %       Temperature values (°C) for the isothermals to draw
    %   'isochorics':    numeric vector
    %       Specific volume values (m3/kg) for the isochorics to draw
    %   'TitleText':    character vector
    %       A chart title
    %
    %   Author: Federico Miretti

    properties
        % Public properties
        % Chart axes limits
        sLim = [0 10]
        hLim = [0 4000]
        % Chart major and minor grid steps
        sMajorGridStep = 0.5
        hMajorGridStep = 25
        sMinorGridStep = 0.1
        hMinorGridStep = 5
        % Chart title
        TitleText = []
        % Values for the iso-lines
        isobarics = 10.^(-2:3)
        isothermals = 100:100:800
        isochorics = [0.05 0.1 0.2 0.5 1 3 8 20]
        isoquality = []
        % List of markers
        markers = struct('Label', {}, 'h', {}, 's', {})
    end
    properties(Access = private)
        % Line objects for the iso-anything lines
        isobaric_lines = {};
        isothermal_lines = {};
        isochoric_lines = {};
        isoquality_lines = {};
        % Annotation objects for the iso-anything lines' labels
        isobaric_labels = {};
        isothermal_labels = {};
        isochoric_labels = {};
        isoquality_labels = {};
        % Scatter object for the markers
        markers_dots = [];
        % Annotation object for the markers' labels
        markers_labels = [];
    end
    properties (Constant)
        T_crit = 373.9460; % C
        p_crit = 220.64; % bar
    end

    % Constructor
    methods
        function obj = SteamChart(varargin)

            % Check if XSteam is installed
            try
                XSteam('h_pT', 100, 500);
            catch ME
                switch ME.identifier
                    case 'MATLAB:UndefinedFunction'
                        link = "<a href=""https://www.mathworks.com/matlabcentral/fileexchange/9817-x-steam-thermodynamic-properties-of-water-and-steam"">install it</a>";
                        error("XSteam cannot be found in the current path. Make sure to " + link + " and to add it to the current path.")
                end
            end

            % Call superclass constructor method
            obj@matlab.graphics.chartcontainer.ChartContainer(varargin{:});
        end
    end

    methods(Access = protected)

        function setup(obj)
            % The setup method runs once when an object is constructed
            
            % Draw saturation lines
            % Temperature breakpoints
            T0 = 0;
            T1 = obj.T_crit;
            x = 1./(1+exp(-10*linspace(-1,1,101)));
            T = T0 + (T1-T0).*x;
            % Evaluate h,s for the saturation lines
            h_sat_vap = arrayfun(@(y) XSteam('hV_T', y), T);
            s_sat_vap = arrayfun(@(y) XSteam('sV_T', y), T);
            h_crit = XSteam('h_pT', obj.p_crit, obj.T_crit);
            s_crit = XSteam('s_pT', obj.p_crit, obj.T_crit);
            h_sat_liq = arrayfun(@(y) XSteam('hL_T', y), T);
            s_sat_liq = arrayfun(@(y) XSteam('sL_T', y), T);
            
            % Get the axes
            ax = getAxes(obj);
            hold(ax,'on')
            
            % Plot saturation lines and crit point
            plot(ax, [s_sat_liq s_sat_vap(end:-1:1)] , [h_sat_liq h_sat_vap(end:-1:1)], 'k', 'LineWidth', 1);
            plot(ax, s_crit, h_crit, '.k', 'MarkerSize', 10)
            
            % Set axis labels
            xlabel(ax, 'Specific Entropy, s (kJ/kg)')
            ylabel(ax, 'Specific Enthalpy, h (kJ/kg)')
            
            % Show grid
            grid(ax, 'on')
            grid(ax, 'minor')
            % Minor grid settings
            ax.MinorGridLineStyle = '-';
            ax.MinorGridAlpha = 1;
            ax.MinorGridColor = [0.9 0.9 0.9];
            % Major grid settings
            ax.GridAlpha = 1;
            ax.GridColor = [0.8 0.8 0.8];
            
            hold(ax,'off')
        end
        
        function update(obj)
            % The update method runs any time the properties are modified
            
            ax = getAxes(obj);
            hold(ax,'on')
            
            % Axis limits
            xlim(ax, obj.sLim)
            ylim(ax, obj.hLim)
            
            % Update grid
            ax.XAxis.MinorTickValues = ax.XAxis.Limits(1):obj.sMinorGridStep:ax.XAxis.Limits(2);
            ax.YAxis.MinorTickValues = ax.YAxis.Limits(1):obj.hMinorGridStep:ax.YAxis.Limits(2);
            ax.XAxis.TickValues = ax.XAxis.Limits(1):obj.sMajorGridStep:ax.XAxis.Limits(2);
            ax.YAxis.TickValues = ax.YAxis.Limits(1):obj.hMajorGridStep:ax.YAxis.Limits(2);
            
            % Update title
            title(getAxes(obj), obj.TitleText);
            
            % Redraw
            obj.draw_isobarics(obj.isobarics);
            obj.draw_isothermals(obj.isothermals);
            obj.draw_isochorics(obj.isochorics);
            obj.draw_isoquality(obj.isoquality);

            obj.draw_markers();
            
            hold(ax,'off')
        end

        function draw_isobarics(obj, p)
            % draw_isobarics(pressure) draws isobarics for the pressure
            % values contained in the pressure vector
            
            % Wipe current isobarics
            for n = 1:length(obj.isobaric_lines)
                delete(obj.isobaric_lines{n});
                delete(obj.isobaric_labels{n});
            end
            
            % get axes object
            ax = getAxes(obj);
            % Curve for labels placement
            ecc = 2;
            x_lim = obj.sLim;
            y_lim = obj.hLim;
            Rx = ecc/(x_lim(2) - x_lim(1));
            Ry = ecc/(y_lim(2) - y_lim(1));
            hyperbole = @(x) (1./((x - x_lim(2)).*Rx.*Ry) + y_lim(2)) .* (x < x_lim(2));
            
            % Draw lines
            s = 0:.1:10;
            for k = 1:length(p)
                h = arrayfun(@(y) XSteam('h_ps', p(k), y), s);
                % Draw line
                obj.isobaric_lines{k} = plot(ax, s, h, 'Color', [0 0 1]);
                % Labels
                [s_label, h_label, orient] = obj.labels_placement(hyperbole, s, h);
                obj.isobaric_labels{k} = text(ax, s_label, h_label, ...
                    ['p = ' num2str(p(k)) ' bar'], 'Color', [0 0 1], ...
                    'Rotation', orient,  'BackgroundColor', 'w', 'Margin', 1);
            end
        end
        
        function draw_isothermals(obj, T)
            % draw_isothermals(temperature) draws isothermals for the
            % temperature values contained in the temperature vector
            
            % Wipe current isothermals
            for n = 1:length(obj.isothermal_lines)
                delete(obj.isothermal_lines{n});
                try
                    delete(obj.isothermal_labels{n});
                end
            end
            
            % get axes object
            ax = getAxes(obj);
            % Curve for labels placement
            ecc = 4.4;
            x_lim = obj.sLim;
            y_lim = obj.hLim;
            Rx = ecc/(x_lim(2) - x_lim(1));
            Ry = ecc/(y_lim(2) - y_lim(1));
            hyperbole = @(x) (1./((x - x_lim(2)).*Rx.*Ry) + y_lim(2)) .* (x < x_lim(2));
            
            p = 10.^(-2:0.01:3);
            for k = 1:length(T)
                h = arrayfun(@(y) XSteam('h_pT', y, T(k)), p);
                s = arrayfun(@(y) XSteam('s_pT', y, T(k)), p);
                % Remove NaNs
                nansIdx = isnan(s) | isnan(h);
                h(nansIdx) = [];
                s(nansIdx) = [];
                obj.isothermal_lines{k} = plot(ax, s, h, 'Color', [220/255,20/255,60/255]);
                % Labels
                s_label_breakpoints = linspace(x_lim(1), x_lim(2), length(s));
                h = interp1(s, h, s_label_breakpoints);
                % Remove NaNs
                nansIdx = isnan(s) | isnan(h);
                h(nansIdx) = [];
                s_label_breakpoints(nansIdx) = [];
                % Place labels
                [s_label, h_label, orient] = obj.labels_placement(hyperbole, s_label_breakpoints, h);
                obj.isothermal_labels{k} = text(ax, s_label, h_label, ...
                    ['t = ' num2str(T(k)) ' °C'], ...
                    'Color', [220/255,20/255,60/255], 'Rotation', orient,  ...
                    'BackgroundColor', 'w', 'Margin', 1);
            end
        end
        
        function draw_isochorics(obj, v)
            % draw_isochore(volume) draws isoisochores for the specific
            % volume values contained in the volume vector
            
            % Wipe current isochorics
            for n = 1:length(obj.isochoric_lines)
                delete(obj.isochoric_lines{n});
                delete(obj.isochoric_labels{n});
            end
            
            % get axes object
            ax = getAxes(obj);
            % Curve for labels placement
            ecc = 5;
            x_lim = obj.sLim;
            y_lim = obj.hLim;
            Rx = ecc/(x_lim(2) - x_lim(1));
            Ry = ecc/(y_lim(2) - y_lim(1));
            hyperbole = @(x) (1./((x - x_lim(2)).*Rx.*Ry) + y_lim(2)) .* (x < x_lim(2));
            
            rho = 1./ v;
            p = [10.^(-2:0.001:1) 10.^(1:0.01:3)];
            for k = 1:length(rho)
                h = arrayfun(@(y) XSteam('h_prho', y, rho(k)), p);
                s = arrayfun(@(y,z) XSteam('s_ph', y, z), p, h);
                % Draw line
                obj.isochoric_lines{k} = plot(ax, s, h, 'Color', [50/255,205/255,50/255]);
                % Draw label
                [s_label, h_label, orient] = obj.labels_placement(hyperbole, s, h);
                obj.isochoric_labels{k} = text(ax, s_label, h_label, ...
                    ['v = ' num2str(v(k)) ' m^3/kg'], ...
                    'Color', [50/255,205/255,50/255], 'Rotation', orient,  ...
                    'BackgroundColor', 'w', 'Margin', 1);
            end
        end
        
        function draw_isoquality(obj, q)
            % draw_isoquality(vapor_quality) draws constant vapor quality
            % lines for the specified values

            % Wipe current isoquality
            for n = 1:length(obj.isoquality_lines)
                delete(obj.isoquality_lines{n});
                delete(obj.isoquality_labels{n});
            end

            % get axes object
            ax = getAxes(obj);
            % Curve for labels placement
            a = 0.3;
            b = 0.3;
            x_lim = obj.sLim;
            y_lim = obj.hLim;
            a = a * ( x_lim(2) - x_lim(1) );
            b = b * ( y_lim(2) - y_lim(1) );
            A = a^2;
            B = - y_lim(1) * a^2;
            C = @(x) ( y_lim(1) * a )^2 + ( x - x_lim(1) ).^2 * b^2 - a^2 * b^2;
            ellipse = @(x) ( (- B + sqrt( B^2 - A*C(x) ) ) ./ A ) ...
                .* ( x <= (x_lim(1)+a) ) .* ( x >= x_lim(1) ) ...
                + ( y_lim(1) + b ) .* ( x < x_lim(1) );

            % Draw lines
            p = 10.^(-2:0.01:3);
            for k = 1:length(q)
                h = arrayfun(@(y) XSteam('h_px', y, q(k)), p);
                s = arrayfun(@(y,z) XSteam('s_ph', y, z), p, h);
                % Remove NaNs
                nansIdx = isnan(s) | isnan(h);
                h(nansIdx) = [];
                s(nansIdx) = [];
                % Draw line
                obj.isoquality_lines{k} = plot(ax, s, h, 'Color', [0.4 0.4 0.4]);
                % Labels
                [s_label, h_label, orient] = obj.labels_placement(ellipse, s, h);
                obj.isoquality_labels{k} = text(ax, s_label, h_label, ...
                    ['x = ' num2str(q(k))], 'Color', [0.4 0.4 0.4], ...
                    'Rotation', orient, 'BackgroundColor', 'w', 'Margin', 1);
            end

        end

        function draw_markers(obj)
            % Draw all markers
            
            ax = getAxes(obj);
            
            if ~isempty(obj.markers_labels)
                delete(obj.markers_dots);
                delete(obj.markers_labels);
            end

            obj.markers_dots = scatter(ax, [obj.markers(:).s], [obj.markers(:).h], 100, 'r', '.');
            obj.markers_labels = text(ax, [obj.markers(:).s]-0.05, [obj.markers(:).h]+20,  {obj.markers(:).Label}, ...
                'Color', 'red', 'FontSize', 14, 'HorizontalAlignment', 'center');

        end
        
        function [s_label, h_label, orient] = labels_placement(obj, labels_guide, s, h)
            % labels_guide is an anonymous function which defines the guide
            % line for the labels. Labels are placed at the intersection
            % between the guide line and the iso-something line they refer to.
            % s and h define the iso-something line
            
            % Find closest intersection
            F = labels_guide(s) - h;
            [~, idx] = min(abs(F));
            % Position and orient the label
            if idx>1 && idx<length(F)
                s_label = s(idx-1) + (s(idx) - s(idx-1)) * (-F(idx-1)) / (F(idx)-F(idx-1));
                h_label = labels_guide(s_label);
                % Orientation
                orient = rad2deg(atan( ( h(idx+1) - h(idx-1) ) / ( s(idx+1) - s(idx-1) ) / diff(obj.hLim) * diff(obj.sLim)));
                % Validate intersection
                out_of_bounds = h_label < obj.hLim(1) || h_label > obj.hLim(2) || s_label < obj.sLim(1) || s_label > obj.sLim(2);
                inaccurate = abs(F(idx)) > ( obj.hLim(2) - obj.hLim(1) ) .* 0.05;
                if out_of_bounds || inaccurate
                    s_label = [];
                    h_label = [];
                    orient = [];
                end
            else
                s_label = [];
                h_label = [];
                orient = [];
            end
        end
        
    end
    
    methods
        function add_marker(obj, label, mode, prop)
            arguments
                obj
                label string {mustBeText}
                mode char {mustBeTextScalar, mustBeMember(mode, {'hs','pT','ph','ps','Lp','LT','Vp','VT'})} = 'hs'
            end
            arguments (Repeating)
                prop
            end
            % add_markers Add markers to the markers list
            %   add_marker(label, 'hs', h, s) adds a labeled marker
            %       specifiying enthalpy (kJ/kg) and enthropy (kJ/kg) coords
            %   add_marker(label, 'pT', p, T) adds a labeled marker
            %       specifiying pressure (bar) and temperature (°C) coords
            %   add_marker(label, 'ph', p, h) adds a labeled marker
            %       specifiying pressure (bar) and enthalpy (kJ/kg) coords
            %   add_marker(label, 'Lp', p) adds a labeled marker on the
            %       lower saturation curve specifiying pressure (bar)
            %   add_marker(label, 'LT', T) adds a labeled marker on the
            %       lower saturation curve specifiying temperature (°C)
            %   add_marker(label, 'Vp', p) adds a labeled marker on the
            %       upper saturation curve specifiying pressure (bar)
            %   add_marker(label, 'VT', T) adds a labeled marker on the
            %       upper saturation curve specifiying temperature (°C)
            %
            %   To add several markers at once, specify label as cell array of
            %   character vectors and physical quantities as vectors

            num_markers = length(label);
            
            switch mode
                case 'hs'
                    h = prop{1};
                    s = prop{2};
                case 'pT'
                    p = prop{1};
                    T = prop{2};
                    s = arrayfun(@(x,y) XSteam('s_pT', x, y), p, T);
                    h = arrayfun(@(x,y) XSteam('h_pT', x, y), p, T);
                case 'ph'
                    p = prop{1};
                    h = prop{2};
                    s = arrayfun(@(x,y) XSteam('s_ph', x, y), p, h);
                case 'ps'
                    p = prop{1};
                    s = prop{2};
                    h = arrayfun(@(x,y) XSteam('h_ps', x, y), p, s);
                case 'Lp'
                    p = prop{1};
                    h = arrayfun(@(y) XSteam('hL_p', y), p);
                    s = arrayfun(@(y) XSteam('sL_p', y), p);
                case 'LT'
                    T = prop{1};
                    h = arrayfun(@(y) XSteam('hL_T', y), T);
                    s = arrayfun(@(y) XSteam('sL_T', y), T);
                case 'Vp'
                    p = prop{1};
                    h = arrayfun(@(y) XSteam('hV_p', y), p);
                    s = arrayfun(@(y) XSteam('sV_p', y), p);
                case 'VT'
                    T = prop{1};
                    h = arrayfun(@(y) XSteam('hV_T', y), T);
                    s = arrayfun(@(y) XSteam('sV_T', y), T);
            end
            
            for n = 1:num_markers
                new_marker.Label = label{n};
                new_marker.h = h(n);
                new_marker.s = s(n);
                obj.markers(end+1) = new_marker;
            end
        end
        
        function remove_marker(obj, label)
            arguments
                obj
                label {mustBeText}
            end
            % remove_markers Remove markers from the markers list
            %   remove_markers(label) Removes the marker with the specified
            %   label. To remove several markers at once, specify label as 
            %   cell array of character vectors
            
            if ~iscell(label)
                label = {label};
            end
            for n = 1:length(label)
                del = ismember({obj.markers(:).Label} , label{n});
                obj.markers(del) = [];
            end
            
            
        end
        
        function save_chart(obj, filename, opts)
            arguments
                obj
                filename {mustBeTextScalar}
                opts.format {mustBeTextScalar} = 'a4'
                opts.orientation {mustBeTextScalar, mustBeMember(opts.orientation, {'portrait','landscape'})} = 'portrait'
            end
            %save_chart Export chart to pdf
            %   save_chart(filename) exports a the chart in a pdf named by
            %   the string or character vector filename
            
            update(obj)
            ax = getAxes(obj);

            fig_export = figure('Visible', 'off');
            copyobj(ax, fig_export);
            
            ax_export = fig_export.Children;
            set(ax_export.Children, 'Clipping', 'on')
            
            fig_export.PaperType = opts.format;
            fig_export.PaperOrientation = opts.orientation;
            
            fig_export.PaperUnits = 'normalized';
            fig_export.PaperPosition = [-0.05 -0.07 1.1 1.14];
            
            saveas(fig_export, filename)
            
        end
        
        function varargout = xlim(obj, varargin)
            if length(varargin) == 1
                obj.sLim = varargin{1};
            elseif isempty(varargin)
                varargout{1} = obj.sLim;
            end
        end
        
        function varargout = ylim(obj, varargin)
            if length(varargin) == 1
                obj.hLim = varargin{1};
            elseif isempty(varargin)
                varargout{1} = obj.hLim;
            end
        end
        
        function title(obj,txt)
            obj.TitleText = txt;
        end
    end
end