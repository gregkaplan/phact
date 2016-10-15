function C=prodt(A,B,nA,nB,ndA)
%function C=prodt(A,B,nA,nB,nsubA)
% Tensor product of A and B along dimensions nA, nB
% argument nsubA is optional.  It is the number of subscripts
% in A.  Without it, if there are trailing singleton dimensions in A beyond nA,
% they are lost from C.
dA=size(A);
if nargin<5
   ndA=ndims(A);
   if nA>ndA
      dA=[dA ones(1,nA-ndA)];
   end
else
   if ndA>ndims(A)
      dA=[dA ones(1,ndA-ndims(A))];
   end
end
dB=size(B);
ndB=ndims(B);
if nB>ndB
   dB=[dB ones(1,nB-ndB)];
end
if ~(dA(nA)==dB(nB)), error('dimension mismatch');end
A=permute(A,[1:(nA-1),(nA+1):ndA nA]);
A=reshape(A,prod(dA([1:nA-1 nA+1:ndA])),dA(nA));
B=permute(B,[nB,1:(nB-1),(nB+1):ndB]);
dA(nA)=[];
dB(nB)=[];
C=zeros([dA dB]);
C(:)=A*B(:,:);
