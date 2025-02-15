function varargout = conic_option_type(varargin)
    %CONIC_OPTION_TYPE [INTERNAL] 
    %
    %  char = CONIC_OPTION_TYPE(char name, char op)
    %
    %Get type info for a particular option.
    %
    %Extra doc: https://github.com/casadi/casadi/wiki/L_1el
    %
    %Doc source: 
    %https://github.com/casadi/casadi/blob/develop/casadi/core/conic.hpp#L560
    %
    %Implementation: 
    %https://github.com/casadi/casadi/blob/develop/casadi/core/conic.cpp#L560-L562
    %
    %
    %
  [varargout{1:nargout}] = casadiMEX(823, varargin{:});
end
