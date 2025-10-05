function varargout = plotcmap(varargin)
%PLOTCMAP plots a colormapped 2-D or 3-D line that maps the values of
% another variable onto the colormap on the line.
% Specifies LineStyle, Marker, and LineWidth in the same way as plot/plot3.
% -------------------------------------------------------------------------
% 
%Syntax:
%
% plotcmap(X, Y, cmap, LineValue)
% plotcmap(Y, cmap)
% plotcmap(X, Y, Z, cmap)
% plotcmap(__, LineValue)
% plotcmap(__, LineSpec)
% plotcmap(__, 'MatchMarkerFaceColor')
% plotcmap(ax, __)
% h = plotcmap(__)
% [h, cb] = plotcmap(__)
% -------------------------------------------------------------------------
%
% Description:
%
% plotcmap(X, Y, cmap, LineValue): plots a 2-D line using plot(X, Y) and
%                                  creates a colormap based on cmap on the
%                                  line which maps the values of LineValue
%                                  onto the colormap.
%                                  X, Y, and LineValue must be mumeric
%                                  vectors with same length.
%                                  cmap can be one of the MATLAB predefined
%                                  colormap names or an customized N X 3
%                                  colormap array.
%                                  (list of predefined colormap names:
%                                  parula, turbo, hsv, hot, cool, spring,
%                                  summer, autumn, winter, gray, bone,
%                                  copper, pink, jet, lines, colorcube,
%                                  prism, flag, white)
%
% plotcmap(X, Y, cmap): plots a 2-D line using plot(X, Y) and creates a
%                       colormap based on cmap on the line without the
%                       values of LineValue where plotcmap maps an implicit
%                       set of values ranges from 1 to length(X).
%
% plotcmap(__, LineSpec): plots a colormapped line with specified
%                         LineStyle, Marker, and LineWidth. Use of LineSpec
%                         is the same as in plot or plot3.
%                         e.g., plotcmap(X, Y, cmap, LineValue, '--o')
%                         plotcmap(X, Y, cmap, 'LineWidth', 2).
%                         Note that plotcmap overwrites the line color
%                         property even if it is specified.
%
% plotcmap(__, 'MatchMarkerFaceColor'): matchs marker face colors of
%                                       each marker to their corresponding
%                                       LineValue to create filled markers.
%
% plotcmap(Y, cmap, __): plots a 2-D colormapped line using plot(Y).
%
% plotcmap(X, Y, Z, cmap, __): plots a 3-D colormapped line using
%                              plot3(X, Y, Z).
%
% plotcmap(ax, __): displays the plot in the target axes. Specify the axes
%                   as the first argument in any of the previous syntaxes.
%                   Same as plot or plot3.
%
% h = plotcmap(__): returns an array of Line objects.
%
% [h, cb] = plotcmap(__): h returns an array of Line objects. Displays
%                         the colorbar based on
%                         caxis([min(LineValue), max(LineValue)]), and cb
%                         returns the colorbar handle.
% -------------------------------------------------------------------------
%
% Examples:
%
% plotcmap(X, Y, parula, LineValue, '--', 'LineWidth', 2):
%   plots a 2-D colormapped line that maps the values of vector LineValue
%   onto the predefined colormap parula with dashed line style and
%   linewidth of 2pt.
%
% plotcmap(X, Y, jet(20), 'o', 'MatchMarkerFaceColor'):
%   plots a 2-D colormapped line that maps an implicit set of values ranges
%   from 1 to length(X) onto 20 colors of the predefined colormap jet with
%   filled circle marker.
%
% plotcmap(ax1, Y, cmap, LineValue):
%   plots a 2-D colormapped line on axis ax1 based on an implicit set of
%   x-coordinate values ranges from 1 to length(Y) that maps the values of
%   vector LineValue onto the customized N X 3 colormap array cmap.
%
% h = plotcmap(X, Y, Z, cool, LineValue, ':s', 'MatchMarkerFaceColor'):
%   plots a 3-D colormapped line that maps the values of vector LineValue
%   onto the predefined colormap cool with dotted line style and filled
%   suqare marker, and returns the line object array to h.
%
% [h, cb] = plotcmap(X, Y, summer, LineValue):
%   plots a 2-D colormapped line that maps the values of vector LineValue
%   onto the predefined colormap summer, returns the line object array to
%   h, displays the colorbar, and returns the colorbar handle to cb.
%
% [~, ~] = plotcmap(X, Y, summer, LineValue):
%   plots a 2-D colormapped line that maps the values of vector LineValue
%   onto the predefined colormap summer and displays the colorbar.
% =========================================================================
%
% version 1.2.0
%   - Added an output argument that displays the colorbar and returns the
%     colorbar handle. E.g., [~, cb] = plotcmap(__).
% Xiaowei he
% 05/22/2022
% -------------------------------------------------------------------------
% version 1.1.1
%   - Fixed a bug that causes errors when plotting on target axes.
%   - Minor updates in headline description.
% Xiaowei He
% 05/07/2022
% -------------------------------------------------------------------------
% version 1.1.0
%   - Added support for using MATLAB predefined colormaps.
%   - Updates in the headline description including a few examples.
% Xiaowei He
% 05/01/2022
% -------------------------------------------------------------------------
% version 1.0.0
% Xiaowei He
% 04/24/2022
% =========================================================================
    % input and output argument check
    nargoutchk(0, 2)
    if nargin < 2
        error('plotcmap(X, Y, Z, cmap, LineValue, LineSpecs): not enough input arguments, must specify at least Y vector and cmap array.')
    else
        % axis handle switch
        axflag = false;
        if ishandle(varargin{1})
            if nargin < 3
                error('plotcmap(ax, X, Y, Z, cmap, LineValue, LineSpecs): not enough input arguments, must specify at least Y vector and cmap array.')
            end
            axflag = true;
        end
    end
    % assign input arguments
    if ~axflag % without axis handle
        iax = 0;
        axmsg = '';
    else % with axis handle
        iax = 1;
        axmsg = 'ax, ';
    end
    % assign cmap
    for i = 1 + iax : nargin
        if ~isvector(varargin{i})
            iarg = i;
            break
        else
            iarg = 0;
        end
    end
    if iarg == 0
        error(['plotcmap(' axmsg '__, cmap): cmap must be an N X 3 numeric array.'])
    else
        cmap = varargin{iarg};
    end
    % check cmap format
    if ~(isnumeric(cmap) && ismatrix(cmap) && size(cmap, 2) == 3)
        error(['plotcmap(' axmsg '__, cmap): cmap must be an N X 3 numeric array.'])
    end    
    % assign X, Y, and Z
    switch iarg
        case 2 + iax % Y
            dflag = 2;
            y = varargin{1+iax};
            if isnumeric(y) && isvector(y)
                x = (1:length(y))';
            else
                error(['plotcmap(' axmsg 'Y, cmap): Y must be a numeric vector.'])
            end
        case 3 + iax % X, Y
            dflag = 2;
            x = varargin{1+iax};
            y = varargin{2+iax};
            if ~(isnumeric(x) && isvector(x) ...
                 && isnumeric(y) && isvector(y) ...
                 && length(x) == length(y))
                error(['plotcmap(' axmsg 'X, Y, cmap): X and Y must be numeric vectors and have the same length.'])
            end
        case 4 + iax % X, Y, Z
            dflag = 3;
            x = varargin{1+iax};
            y = varargin{2+iax};
            z = varargin{3+iax};
            if ~(isnumeric(x) && isvector(x) ...
                 && isnumeric(y) && isvector(y) ...
                 && isnumeric(z) && isvector(z) ...
                 && length(x) == length(y) && length(x) == length(z))
                error(['plotcmap(' axmsg 'X, Y, Z, cmap): X, Y, and Z must be numeric vectors and have the same length.'])
            end
        otherwise
            error(['plotcmap(' axmsg 'X, Y, Z, cmap): too many input arguments.'])
    end
    % LineValue and line properties
    argflag = false;
    mmfcflag = false;
    if length(varargin) > iarg 
        argflag = true;
        if strcmp(varargin{end}, 'MatchMarkerFaceColor') % MatchMarkerFaceColor switch
            mmfcflag = true;
            argtemp = varargin(iarg+1:end-1);
        else
            argtemp = varargin(iarg+1:end);
        end
    end  
    if argflag
        if isnumeric(argtemp{1}) % LineValue switch on
            v = argtemp{1};
            if ~(isvector(v) && length(v) == length(y)) % check LineValue format
                error('plot(__, LineValue): LineValue must a vector and have the same length as Y.')
            end
            if length(argtemp) > 1 % check if there are line property arguments after LineValue
                LineSpec = argtemp(2:end);
            else
                LineSpec{1} = char.empty;
            end
        else % LineValue switch off
            v = (0 : length(x)-1);
            LineSpec = argtemp;
        end
    else % no LineValue or line property
        v = (0 : length(x)-1);
        LineSpec{1} = char.empty;
    end
    % map LineValue to cmap
    vfrac = (v - min(v))/(max(v) - min(v)); % v(i) location in v vector
    ic = round(vfrac*(size(cmap, 1) - 1)) + 1; % index in cmap corresonding to v(i)
    % plot
    % determine plot axis
    if axflag
        ax = varargin{1};
    else
        ax = gca;
    end
    holdflag = ishold; % current plot hold status
    hold on
    L = length(x);
    h = gobjects(L-1, 1); % initiate line handle array
    switch dflag
        case 2 % 2-D plot
            for i = 1 : L-1
                h(i) = plot(ax, [x(i), x(i+1)], [y(i), y(i+1)], ...
                            LineSpec{:}, 'Color', cmap(ic(i), :));
                if mmfcflag
                    set(h(i), 'MarkerFaceColor', cmap(ic(i), :))
                end
            end
        case 3 % 3-D plot
            for i = 1 : L-1
                h(i) = plot3(ax, [x(i), x(i+1)], [y(i), y(i+1)], [z(i), z(i+1)], ...
                             LineSpec{:}, 'Color', cmap(ic(i), :));
                if mmfcflag
                    set(h(i), 'MarkerFaceColor', cmap(ic(i), :))
                end
            end
            view(-15, 20)
    end
    % colormap corection if plotting markers only
    if strcmp(h(1).LineStyle, 'none') && ~strcmp(h(1).Marker, 'none')
        delete(h(:))
        h = gobjects(L, 1);
        switch dflag
            case 2 % 2-D plot
                for i = 1 : L
                    h(i) = plot(ax, x(i), y(i), ...
                                LineSpec{:}, 'Color', cmap(ic(i), :));
                    if mmfcflag
                        set(h(i), 'MarkerFaceColor', cmap(ic(i), :))
                    end
                end
            case 3 % 3-D plot
                for i = 1 : L
                    h(i) = plot3(ax, x(i), y(i), z(i), ...
                                LineSpec{:}, 'Color', cmap(ic(i), :));
                    if mmfcflag
                        set(h(i), 'MarkerFaceColor', cmap(ic(i), :))
                    end
                end                
        end
    end
    % restore plot hold status
    if ~holdflag
        hold off
    end
    % output
    if nargout == 1
        varargout{1} = h;
    elseif nargout == 2
        colormap(ax, cmap)
        caxis([min(v), max(v)])
        cb = colorbar(ax);
        varargout{1} = h;
        varargout{2} = cb;
    end
endf