function [A,B,Q,Z,v] = qzdiv(stake,A,B,Q,Z,v)
%function [A,B,Q,Z,v] = qzdiv(stake,A,B,Q,Z,v)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of abs(B(i,i)/A(i,i))>stake are in lower right 
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.  The columns of v are sorted correspondingly.
%
% by Christopher A. Sims
% modified (to add v to input and output) 7/27/00
vin = nargin==6;
if ~vin, v=[]; end;
[n jnk] = size(A);
root = abs([diag(A) diag(B)]);
root(:,1) = root(:,1)-(root(:,1)<1.e-13).*(root(:,1)+root(:,2));
root(:,2) = root(:,2)./root(:,1);
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if (root(j,2) > stake | root(j,2) < -.1) 
         m=j;
         break
      end
   end
   if (m==0) 
      return 
   end
   for k=m:1:i-1
      [A B Q Z] = qzswitch(k,A,B,Q,Z);
      tmp = root(k,2);
      root(k,2) = root(k+1,2);
      root(k+1,2) = tmp;
      if vin
         tmp=v(:,k);
         v(:,k)=v(:,k+1);
         v(:,k+1)=tmp;
      end
   end
end         
