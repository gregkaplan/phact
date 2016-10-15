function varargout = size(x, varargin)
% by SeHyoun Ahn, July 2016
varargout = {size(x.values,varargin{:})};
