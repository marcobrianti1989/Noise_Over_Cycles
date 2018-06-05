%-----------------------------------------------------------------------------
% Compute Impulse Responses and A0 matrices for 
% the penalty function approach under a general fiscal rule.
%-----------------------------------------------------------------------------

function [ffactor, A0PF, Ltilde] = penalty_fun_general(Sigma,n,p,F,J,Hor,nperUhlig)

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
    nper = nperUhlig;
    IR = LtildeTemp;

    %------------------------------
    % Identify business cycle shock          
    %------------------------------
    
    Aeq=[];
    beq=[];
    q1ga=rand(n,1);
    [Qhat1,~] = fmincon(@penalty_bc_shock,q1ga,[],[],Aeq,beq,[],[],@mycon,optimset('MaxFunEvals',40000,'MaxIter',20000,'Display','off','Algorithm','active-set'));
    
    %-------------------------------
    % Identify monetary policy shock          
    %-------------------------------
    
    Aeq=Qhat1';
    beq=0;
    q2ga=rand(n,1);
    [Qhat2,~] = fmincon(@penalty_mp_shock,q2ga,[],[],Aeq,beq,[],[],@mycon,optimset('MaxFunEvals',40000,'MaxIter',20000,'Display','off','Algorithm','active-set'));
    
    %------------------------------
    % Identify inflation shock          
    %------------------------------
    
    Z = zeros(n-3,n);
    Z(:,1:2) = eye(2);
    Rx = [Z*LC; Qhat1'; Qhat2'];
    Nx=null(Rx);
    Xhat=randn(n,1);
    Qhat3 =Nx*Nx'*Xhat/norm(Nx'*Xhat);
    
    %-----------------------------------------------------------------------------------------------------
    % Identify spending shock (under assumption that spending does not respond to contemporaneously taxes)
    %-----------------------------------------------------------------------------------------------------
    
    Z = zeros(n-4,n);
    Z(:,1) = 1;
    Rx = [Z*A0; Qhat1';Qhat2';Qhat3'];
    Nx=null(Rx);
    Xhat=randn(n,1);
    Qhat4 =Nx*Nx'*Xhat/norm(Nx'*Xhat);
    
    %-----------------------------------------------------------------------------------------------------
    % Identify tax shock (uniquely identified by orthogonality to the other four shocks)
    %-----------------------------------------------------------------------------------------------------
    
    Rx = [Qhat1';Qhat2';Qhat3';Qhat4'];
    Nx=null(Rx);
    Xhat=randn(n,1);
    Qhat5 =Nx*Nx'*Xhat/norm(Nx'*Xhat);
    
    %---------------------------------------------
    % Calculate elasticities and impulse responses
    %---------------------------------------------

    Q = [Qhat1 Qhat2 Qhat3 Qhat4 Qhat5]; 
    ffactor = LC*Q;
    A0PF    = A0*Q;

    Omega1 = [ffactor;zeros((p-1)*n,size(ffactor,2))];
    Ltilde  = vm_irf(F,J,ffactor,Hor,n,Omega1);

end
