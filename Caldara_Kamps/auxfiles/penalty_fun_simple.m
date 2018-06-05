%-----------------------------------------------------------------------------
% Compute Impulse Responses and A0 matrices for 
% the penalty function approach under a simple fiscal rule.
%-----------------------------------------------------------------------------

function [ffactor, A0PF, Ltilde] = penalty_fun_simple(Sigma,n,p,F,J,Hor,nperPF)

    %--------------
    % Preliminaries
    %--------------
    
    global IR;
    global ssigma;
    global nper;
    LC =chol(Sigma,'lower');
    A0 = (LC')\eye(size(LC,1));

    Omega1 = [LC;zeros((p-1)*n,size(LC,2))];  
    LtildeTemp    = vm_irf(F,J,LC,4,n,Omega1);
    s = sqrt(diag(Sigma));
    ssigma = s;
    nper = nperPF;
    IR = LtildeTemp;
    
    %------------------------------
    % Identify business cycle shock          
    %------------------------------
    
    Aeq=[];
    beq=[];
    q1ga=rand(n,1);
    [Qhat1,~] = fmincon(@penalty_bc_shock,q1ga,[],[],Aeq,beq,[],[],@mycon,optimset('MaxFunEvals',40000,'MaxIter',20000,'Display','off','Algorithm','active-set'));

    %---------------------------------------------------------------------------------
    % Identify spending shock (imposing a simple fiscal rule through zero restrictions)
    %---------------------------------------------------------------------------------
    
    Z = zeros(n-2,n); Z(1,1) = 1; Z(2:3,4:5) = eye(2); % Zero restrictions;
    Rx = [Z*A0; Qhat1'];
    Nx=null(Rx);
    Xhat=randn(n,1);
    Qhat2 =Nx*Nx'*Xhat/norm(Nx'*Xhat);
    
    %---------------------------------------------------------------------------------
    % Identify Tax shock  (imposing a simple fiscal rule through zero restrictions)
    %---------------------------------------------------------------------------------
        
    Z = zeros(n-2,n); Z(1,2) = 1; Z(2:3,4:5) = eye(2); % Zero restrictions;
    Rx = [Z*A0; Qhat1'];
    Nx=null(Rx);
    Xhat=randn(n,1);
    Qhat3 =Nx*Nx'*Xhat/norm(Nx'*Xhat);

    %---------------------------------------------
    % Calculate elasticities and impulse responses
    %---------------------------------------------
    
    Q = [Qhat1 Qhat2 Qhat3]; 
    ffactor = LC*Q;
    A0PF    = A0*Q;

    Omega1 = [ffactor;zeros((p-1)*n,size(ffactor,2))];
    Ltilde  = vm_irf(F,J,ffactor,Hor,n,Omega1);

end
