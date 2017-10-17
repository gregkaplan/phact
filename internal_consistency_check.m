function [epsilon] = internal_consistency_check(G1,impact,n_g_red,from_red,to_red,g1,psi,F,n_v,n_g,T,steadystate,plotting,IRF,dt)
% DOCSTRING TO BE UPDATED
% Does an internal consistency check
% XXX A different consistency check will be implemented XXX
%
% by SeHyoun Ahn, March 2017
%
% REFERENCE: Ahn, SeHyoun, Greg Kaplan, Benjamin Moll, Thomas Winberry, and
%    Christian Wolf. "When Inequality Matters for Macro and Macro Matters
%    for Inequality."
%
% PARAMETERS:
%    G1 = dynamics equation of the reduced model (usually output of
%         schur_solver)
%    impact = matrix of shock to variables (usually an output of
%             schur_solver)
%    n_g_red = number of state variables in the reduced model
%    from_red = inverse projection matrix from the reduced basis
%    to_red = projection matrix to the reduced basis
%    g1 = dynamics equation for the full model
%    psi = impact of shocks on state variables
%    F = matrix transform to get stable parts of v (usually output of
%        schur_solver)
%    n_v = number of choice variables in full model
%    n_g = number of state variables in full model
%    T = time period to check consistency with
%    steadystate = steady-state values of each variables
%    plotting = 1 to make diagnostic plots
%    IRF = 1 if compute consistency check for one a period shock
%          0 if compute simulated process
%    dt = (optional) if a good guess of dt is available the step to
%         estimate time-step can be bypassed
%
% OUTPUTS:
%    epsilon = relative errors for internal consistency check
%              (check paper in reference for details)
%
% EXAMPLES:
%     This file requires very specific examples. See the example given in
%     <examples/KrusellSmith/mainfile.m> example (provided from github at
%     <https://github.com/gregpkaplan/phact>)
%
% SYNTAX (you can copy and paste the following) :
% [epsilon] = internal_consistency_check(G1,impact,n_g_red,from_red,to_red,g1,psi,F,n_v,n_g,T,steadystate,plotting,IRF,dt)


n_v_red = size(G1,1) - n_g_red;
n_p = size(g1,1) - n_v - n_g;

% Parse g1 to have easy access
B_pv = -g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,1:n_v);
B_pg = -g1(n_v+n_g+1:n_v+n_g+n_p,n_v+n_g+1:n_v+n_g+n_p)\g1(n_v+n_g+1:n_v+n_g+n_p,n_v+1:n_v+n_g);
B_pvg = [B_pv,B_pg];
B_gg = g1(n_v+1:n_v+n_g,n_v+1:n_v+n_g);
B_gv = g1(n_v+1:n_v+n_g,1:n_v);
B_gp = g1(n_v+1:n_v+n_g,n_v+n_g+1:n_v+n_g+n_p);

% Adjust F with projection matrices
F = @(x) from_red(1:n_v,1:n_v_red)*(F*(to_red(n_v_red+1:n_v_red+n_g_red,n_v+1:n_v+n_g)*x));

% Define Dynamics A given F
A = @(x) B_gg*x + B_gp*(B_pg*x) + B_gv*F(x) + B_gp*(B_pv*F(x));

% Compute time step for explicit update
% eigs can fail to find the largest eigenvalue, so it progressively
% increases the size to find the largest eigenvalue. This part of the code
% might be updated in the future to improve the search speed. (March 2017)
if nargin < 15
    aux = 5;
    while true
        try
            if aux > min(n_g,5000)
                dt = 1e-6;
                warning('eigenvalues not found');
                break;
            end
            dt = 1/max(abs(eigs(A,n_g,aux)))/2;
            break
        catch
            aux = aux+5;
            fprintf('<internal_consistency_check>: Failed to find an approximation for the largest eigenvalue\n');
            fprintf('                              adjusting parameter k of eigs to %d\n',aux);
        end
    end
end
N = ceil(T/dt);
vTime = linspace(0,(N-1)*dt,N);

% Initial distribution for a given impact
if IRF
    shock = zeros(size(impact,2),N+1);
    shock(:,1) = 1;
else
    shock = randn(size(impact,2),N+1);
end
g_red = dt^0.5*impact*shock(:,1);
g = dt^0.5*psi(n_v+1:n_v+n_g)*shock(:,1);

% If step size is too small only keep sampled space
if N > 10000
    step = floor(N/5000);
    p_red = zeros(n_p,size(impact,2),floor((N-1)/step)+1);
    p = zeros(n_p,size(impact,2),floor((N-1)/step)+1);
    vTime = vTime(step*(0:floor((N-1)/step))+1);
else
    step = 1;
    p_red = zeros(n_p,size(impact,2),N);
    p = zeros(n_p,size(impact,2),N);
end
if N>100
    tstart = tic;
end
for i = 1:N
    if mod(i,step) == 0
        p_red(:,:,i/step) = B_pvg*from_red(1:n_v+n_g,:)*g_red+steadystate(n_v+n_g+1:n_v+n_g+n_p);
        p(:,:,i/step) = B_pv*F(g)+B_pg*g+steadystate(n_v+n_g+1:n_v+n_g+n_p);
        %fprintf('%d of %d iterations\n',i,N);
    end
    g_red = g_red + dt*G1*g_red + (dt^(1/2))*impact*shock(:,i+1);
    g = g + dt*A(g) + (dt^(1/2))*psi(n_v+1:n_v+n_g)*shock(:,i+1);
    if i == 100
        timed = toc(tstart);
        estimated_runtime = timed*N/100/60;
        fprintf('<check_internal>: Estimated simulation time is %.6f minutes\n',estimated_runtime);
    end
end

if plotting
    for i = 1:size(p,2)
        aux = reshape((abs(p(:,i,:)-p_red(:,i,:))./abs(p(:,i,:))),size(p,1),size(p,3));
        figure;
        hold on;
        for j = 1:size(p,1)
            plot(vTime,aux(j,:)','DisplayName',['variable ',num2str(j)]);
        end
        title(['Internal Consistency Check for Shock ',num2str(i)]);
        xlabel('time');
        ylabel('relative errors');
        legend('show');
    end
end
epsilon = max(abs(p_red-p)./abs(p),[],3);
fprintf('<internal_consistency_check>: The maximum relative error is %.6e\n',max(max(epsilon)));