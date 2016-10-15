function [A,B,Q,Z] = qzdivct(stake,A,B,Q,Z)
%function [A,B,Q,Z] = qzdivct(stake,A,B,Q,Z)
%
% Takes U.T. matrices A, B, orthonormal matrices Q,Z, rearranges them
% so that all cases of real(B(i,i)/A(i,i))>stake are in lower right
% corner, while preserving U.T. and orthonormal properties and Q'AZ' and
% Q'BZ'.  abs(A(i,i))<1e-11 is interpreted as a zero and as generating
% an infinitely positive real part of the ratio.  All i's for which this
% criterion are satisfied are grouped together in the lower right corner
% of the lower right corner, with the non-zero roots above them.  This
% version differs from
% qzdiv in that it works on the real part's value, as is appropriate for
% continuous time models, instead of on the absolute value, as is
% appropriate for discrete time models.
%
realsmall=sqrt(eps)*10;
%realsmall=1e-3;
[n jnk] = size(A);
root = [diag(A) diag(B)];
% first sort on the non-zero root criterion
xdown0 = abs(root(:,1))<realsmall;
xdown = (xdown0 | (real(root(:,2)./(xdown0+root(:,1))) > stake));
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if xdown0(j)
         m=j;
         break
      end
   end
   if (m==0)
      break
   end
   for k=m:1:i-1
      [A B Q Z] = qzswitch(k,A,B,Q,Z);
		root=[diag(A) diag(B)];
		xdown0(k:k+1)=flipud(xdown0(k:k+1));
		xdown(k:k+1)=flipud(xdown(k:k+1));
		if any(xdown(k:k+1)~=(xdown0(k:k+1) | (real(root(k:k+1,2)./(xdown0(k:k+1)+root(k:k+1,1)))) > stake))
			disp('xdown shift during 0 pack at i,k:')
			disp([i k])
		end
   end
end
% now repeat, using the stake criterion
for i = n:-1:1
   m=0;
   for j=i:-1:1
      if xdown(j)
         m=j;
         break
      end
   end
   if (m==0)
      return
   end
   for k=m:1:i-1
		gevOld=root(k:k+1,:);
		[A B Q Z] = qzswitch(k,A,B,Q,Z);
		root=[diag(A) diag(B)];
		xdown0(k:k+1)=flipud(xdown0(k:k+1));
		xdown(k:k+1)=flipud(xdown(k:k+1));
		if any(xdown(k:k+1)~=(xdown0(k:k+1) | (real(root(k:k+1,2)./(xdown0(k:k+1)+root(k:k+1,1)))) > stake))
			disp('xdown shift during pos pack at i,k:')
			disp([i k])
			gev=root(k:k+1,:);
			[gevOld gevOld(:,1).\gevOld(:,2);gev gev(:,1).\gev(:,2)]
		end
   end
end