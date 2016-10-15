function x = fsolve(h,x0,varargin)
% by SeHyoun Ahn, Jan 2016
% Check included README.pdf to see the syntax for vector valued fsolve
x=myAD(x0);
tmp=@(z) h([z;getvalues(varargin{1})]);
x.values=fsolve(tmp,x0,varargin{2:end});
z=myAD([x.values;getvalues(varargin{1})]);
tmp_deriv=getderivs(h(z));
x.derivatives=-tmp_deriv(:,1:length(x0))\tmp_deriv(:,length(x0)+1:end)*getderivs(varargin{1});
