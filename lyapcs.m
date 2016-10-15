function [F,eno]=lyapcs(A,B,C)
% function [F,eno]=lyapcs(A,B,C)
% solves AF+FB=-C
[ra,ca]=size(A);
[rb,cb]=size(B);
[rc,cc]=size(C);
eno=0;
if ~(ra==ca & rb==cb & rc==ca & cc==rb)
	disp('dimensions conflict')
	eno=1;
else
	[ua,ta]=schur(A);
	[ub,tb]=schur(B);
	[ua,ta]=rsf2csf(ua,ta);
	[ub,tb]=rsf2csf(ub,tb);
	C=ua'*C*ub;
	% matlab's routine just checks for matching diagonal elt's and quits if it
	% finds any.  This routine finds one of the multiple solutions if one exists and
	% returns an F that in some sense tries to solve the equation when an exact solution
	% does not exist
	%---------------------------------------
	%crita=repmat(diag(a),1,rb);
	%critb=repmat(diag(b)',ra,1);
	%if any(any(abs(crita)+abs(critb)<1000*eps))
	%	disp('matching singularities')
	%	eno=2;
	%else
	%-------------------------------------
	F=zeros(ca,rb);
	for cf=1:rb;
		x=ta+tb(1,1)*eye(ra);
		[u,d,v]=svd(x);
		s=sum(diag(d));
		if s<1e7*realmin
			F(:,cf)=zeros(ca,1);
			eno=2;
		else
			dr=diag(d)/s;
			nrank=min(find(dr<eps*1000))-1;
			if isempty(nrank)
				nrank=ca;
			end
			if nrank<ca
				eno=2;
			end
			F(:,cf)=-(v(:,1:nrank).*repmat(1../diag(d(1:nrank,1:nrank))',ca,1))*u(:,1:nrank)'*C(:,1);
		end
		if rb>cf
			C=C(:,2:rb-cf+1)+F(:,cf)*tb(1,2:rb-cf+1);
			tb=tb(2:rb-cf+1,2:rb-cf+1);
		end
	end
	F=ua*F*ub';
	if eno==2
		disp('matched eigenvalues. multiple or nonexistent solutions')
		% it is left to the user whether, and by what criterion, to check AF+FB=C in this case.
	end
end