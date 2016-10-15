function  [GS,pickn,ok,uis,uiu,vs]=gstate(G1,impact,pick)
%function  [GS,pickn,ok,uis,vs]=gstate(G1,impact,pick)
% G1:    the coefficient on  lagged y from the output of gensys.m.  
% impact:the coefficient on exogenous shocks from the output of gensys.m
% pick:  an optional guess at a matrix of coefficients that extracts 
%        the state vector from y
% ok:     0 if pick matrix is just no good.
%         1 if pick matrix is usable forward, but loses initial date information
%         2 if pick matrix is not usable forward, is part of a correct state vector, but is not complete
%         3 if pick matrix is     usable forward, is part of a correct state vector, but is not complete
%         4 if pick matrix is not usable forward, summarizes past, but is redundant
%         5 if pick matrix is     usable forward, summarizes past, but is redundant
%         6 if pick matrix summarizes past and is not redundant, but is not usable forward
%         7 if pick matrix is perfect, both forward and as history summary
% pickn: a matrix of coefficients that extracts the state vector from y.  Equal
%        to pick if pick is supplied and ok=1; otherwise, an ok-maximizing pick matrix that
%        is somewhat similar to pick.             
% GS:    the matrix of coefficients on the lagged state variable
% uis,vs: If uis'*vs is square and full rank, ok=7 is possible, otherwise not.  If ok<2, an ok>2 can be found by trying
%        a pick in the row space of vs'.  Any pick with pick*uis full column rank will provide a foward state 
%        (i.e. ok an odd number).
% uiu:   uiu'y(t)=0, and this is interpretable as the decision rule setting controls as functions of the states.
% The solution was in the form y(t)=G1*y(t-1)+impact*z(t).  Now it's in the form y(t)=GS*pickn*y(t-1)+impact*z(t).
% In general, pickn is mxn with m<n.
REALSMALL=1e-9;
[nr,nc]=size(G1);
%if nr<nc
%   G1=[G1;zeros(nc-nr,nc)];
%end
[u d v]=svd(G1);
top=find(diag(d)>REALSMALL);
nd=top(end);
us=u(:,top);
uu=u(:,top+1:end);
d=d(top,top);
vs=v(:,top);
if nargin<=2
   pick=vs';
   vp=v;
   vps=vs;
   dp=eye(nd);
   ups=eye(nr,nd);
   ndp=nd;
else
   [up dp vp]=svd(pick);
   topp=find(diag(dp)>REALSMALL);
   ndp=topp(end);
   vps=vp(:,topp);
   dp=dp(topp,topp);
   ups=up(:,topp);
end
    % If we were worried about efficiency, we'd skip some of this when pick=vs.
    %Does pick summarize history?  (pick in v' row space, v in vp' row space).
    pinv=all(all((pick'-vs*vs'*pick').^2<REALSMALL));
    vinp=all(all((vs-vps*vps'*vs).^2<REALSMALL));
    okpast=pinv+vinp;
   % Does pick summarize all current info?  (impact in us column space, pick*uu full rank)
   [ui,di,vi]=svd([impact us]);
   topi=find(diag(di)>REALSMALL);
   ndi=length(topi);
   uis=ui(:,topi);
   uiu=ui(:,ndi+1:size(ui,2));
   if ndi<size(G1,1)
      if(size(pick,1)<size(uis,2))
              oknow=0;
      else
          [ut,dt,vt]=svd(pick*uis);
          toppu=find(diag(dt)>REALSMALL);        
          if length(toppu)<size(us,2)
              oknow=0;
          else
              oknow=1;
          end
      end
   else
      oknow=0;
   end
   if vinp
      GS=G1/pick;
      pickn=pick;
   elseif pinv
      r=vs-vps*vps'*vs;
      [ur,dr,vr]=svd(r);
      topr=find(dr>REALSMALL);
      p2=ur(:,topr);
      pickn=[pick;p2'];
      GS=G1/pickn;
   elseif oknow
      GS=G1/[pick;uiu'];
      GS=GS(:,1:size(pick,1));
      pickn=pick;
   else
      pickn=vs';
      GS=us*d;
   end
ok=oknow+2*pinv+4*vinp;
