function varargout = nlpsol_options(varargin)
    %NLPSOL_OPTIONS [INTERNAL] 
    %
    %  {char} = NLPSOL_OPTIONS(char name)
    %
    %Get all options for a plugin.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1t5
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.hpp#L825
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/develop/casadi/core/nlpsol.cpp#L825-L827
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(835, varargin{:});
end
