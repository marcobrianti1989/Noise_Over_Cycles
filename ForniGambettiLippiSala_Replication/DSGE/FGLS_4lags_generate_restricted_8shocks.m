function [CUM_SER,SER,Z_SER,IRF,ZIRF,LIRF] = FGLS_4lags_generate_restricted_8shocks(param,obs,restrict_v);


if restrict_v == 1;
    param(7) = 0;
end

rho_d = param(1);
rho_q = param(2);
rho_p = param(3);
rho_w = param(4);
rho_g = param(5);

sig_eps = param(6);
%sig_eps = (1-rho)*sig_u;
%sig_eta = sig_u*(rho^.5);
sig_v = param(7);
sig_p = param(8);
sig_w = param(9);
sig_g = param(10);
sig_d = param(11);
sig_q = param(12);

h = param(13);
alp = param(14);
phi = param(15);
chi = param(16);
xi = param(17);
cal = param(18);
cal_w = param(19);
iot = param(20);
iot_w = param(21);
mu_p = param(22);
theta_p = param(23); nu_p = theta_p;
mu_w = param(24);
theta_w = param(25); nu_w = theta_w;

rho_r = param(26);
gam_pi = param(27);
gam_y = param(28);

zet = param(29); % from Justiniano, Primiceri, and Tambalotti
del = param(30);
bet = param(31);

sig_ns = param(32);

% shocks: da d q m_p epsma_p m_w epsma_w g
% endogenous: lam psi c i kbar kk u y n rk w pi mc mc_w r ri ca ia ya ka wa dca diad dya dn dwa dr dpi dy;

lam_p=1; psi_p=2; c_p=3; i_p=4; kbar_p=5; kk_p=6; u_p=7; y_p=8; n_p=9; rk_p=10; w_p=11; pi_p=12; mc_p=13; mc_w_p=14; r_p=15;ri_p=16;
n_eq = 16;

FU = zeros(n_eq,n_eq);% FF lead
PR = zeros(n_eq,n_eq);% G  contemp
PA = zeros(n_eq,n_eq);% H  lags

%ca_p=17; ia_p=18;ya_p=19;ka_p=20; wa_p=21; dca_p=22; diad_p=23; dya_p=24; dn_p=25; dwa_p=26; dr_p=27; dpi_p=28; dy_p=29;
%FU = zeros(29,29);% FF lead
%PR = zeros(29,29);% G  contemp
%PA = zeros(29,29);% H  lags

% states

% da = pty - pty(-1)

eps_t = 1;eps1_t = 2;v_t = 3;pty_t = 4;mp_t = 5;e_mp_t = 6;mw_t = 7;e_mw_t = 8;g_t = 9; d_t = 10;q_t = 11;s_t = 12;dpty_t = 13;dpty1_t = 14;eps2_t = 15;eps3_t = 16;
n_states = 16;

FS = zeros(n_eq,n_states);% L
PRS = zeros(n_eq,n_states);% M

% 1st eq:
PR(1,lam_p) = 1;
FU(1,c_p) = - (h*bet/((1-h*bet)*(1-h)));
PR(1,c_p) = ((1+h^2*bet)/((1-h*bet)*(1-h)));
PA(1,c_p) = - (h/((1-h*bet)*(1-h)));

FS(1,dpty_t) = - (h*bet/((1-h*bet)*(1-h)));
PRS(1,dpty_t) = (h/((1-h*bet)*(1-h)));
%lam = (h*bet/((1-h*bet)*(1-h)))*c(+1) - ((1+h^2*bet)/((1-h*bet)*(1-h)))*c + (h/((1-h*bet)*(1-h)))*c(-1) +
% + (h*bet/((1-h*bet)*(1-h)))*da(+1) - (h/((1-h*bet)*(1-h)))*da;

% 2nd eq:
PR(2,lam_p) = 1;
FU(2,lam_p) = -1;
PR(2,ri_p) = -1;

FS(2,dpty_t) = 1;
%lam = lam(+1) - da(+1) + ri;

% 3rd eq:
PR(3,ri_p) = 1;
PR(3,r_p) = -1;
FU(3,pi_p) = 1;
%ri = r - pi(+1);

% 4th eq:
PR(4,psi_p) = 1;
FU(4,psi_p) = -(1-del)*bet;
FU(4,lam_p) = -(1-(1-del)*bet);
FU(4,rk_p) = -(1-(1-del)*bet);

FS(4,dpty_t) = (1-del)*bet + (1-(1-del)*bet) ;
%psi = (1-del)*bet*(psi(+1)-da(+1)) + (1-(1-del)*bet)*(lam(+1)-da(+1)+rk(+1));
%psi = (1-del)*bet*psi(+1) - (1-del)*bet*da(+1) - (1-(1-del)*bet)*da(+1) + (1-(1-del)*bet)*lam(+1) + (1-(1-del)*bet)*rk(+1);

% 5th eq
PR(5,lam_p) = 1;
PR(5,psi_p) = -1;
PR(5,i_p) = chi+bet*chi;
PA(5,i_p) = - chi;
FU(5,i_p) = - bet*chi;

FS(5,dpty_t) = - bet*chi;
PRS(5,dpty_t) = chi;

PRS(5,d_t) = -1;
%lam = psi + d - chi*(i-i(-1)+da) + bet*chi*(i(+1)-i+da(+1));
%lam = psi + d + chi*i(-1) - chi*da + bet*chi*i(+1) - chi*i - bet*chi*i + bet*chi*da(+1);

% 6th eq
PR(6,rk_p) = 1;
PR(6,u_p) = -xi;
%rk = xi*u;

% 7th eq
PR(7,kk_p) = 1;
PR(7,u_p) = -1;
PA(7,kbar_p) = -1;

PRS(7,dpty_t) = 1;
%kk = u + kbar(-1) - da;

% 8th eq
PR(8,kbar_p) = 1;
PA(8,kbar_p) = -(1-del);
PR(8,i_p) = -del;

PRS(8,dpty_t) = (1-del);
PRS(8,d_t) = -del;
% kbar = (1-del)*(kbar(-1)-da) + del*(d+i); kbar = (1-del)*kbar(-1) + del*i - (1-del)*da + del*d

% 9th eq
PR(9,y_p) = 1;
PR(9,kk_p) = -alp;
PR(9,n_p) = -(1-alp);
%y = alp*kk + (1-alp)*n;

% 10th eq
PR(10,pi_p) = 1;
PA(10,pi_p) = -(iot/(1+iot*bet));
FU(10,pi_p) = -(bet/(1+iot*bet));
PR(10,mc_p) = -(((1-cal*bet)*(1-cal))/(cal*(1+iot*bet)));

PRS(10,mp_t) = -1;
% pi = (iot/(1+iot*bet))*pi(-1) + (bet/(1+iot*bet))*pi(+1) + (((1-cal*bet)*(1-cal))/(cal*(1+iot*bet)))*mc + m_p;

% 11th eq
PR(11,mc_p) = 1;
PR(11,rk_p) = -alp;
PR(11,w_p) = -(1-alp);
%mc = alp*rk + (1-alp)*w;

% 12th eq
PR(12,kk_p) = 1;
PR(12,n_p) = -1;
PR(12,w_p) = -1;
PR(12,rk_p) = 1;
%kk - n = w - rk;

% 13th eq
PR(13,w_p) = 1;
PA(13,w_p) = -(1/(1+bet));
FU(13,w_p) = -(bet/(1+bet));

PA(13,pi_p) = -(iot_w/(1+bet));
PR(13,pi_p) = ((1+bet*iot_w)/(1+bet));
FU(13,pi_p) = -(bet/(1+bet));

PR(13,mc_p) = (((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+phi*(1+1/mu_w))));

PRS(13,dpty1_t) = - (iot_w/(1+bet));
PRS(13,dpty_t) = ((1+bet*iot_w)/(1+bet));
FS(13,dpty_t) = - (bet/(1+bet));

PR(13,mc_w_p) = (((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+phi*(1+1/mu_w))));

PRS(13,mw_t) = -1;
%w = (1/(1+bet))*w(-1) + (bet/(1+bet))*w(+1) + (iot_w/(1+bet))*pi(-1) - ((1+bet*iot_w)/(1+bet))*pi + (bet/(1+bet))*pi(+1) +
% (iot_w/(1+bet))*da(-1) - ((1+bet*iot_w)/(1+bet))*da + (bet/(1+bet))*da(+1) - (((1-cal_w*bet)*(1-cal_w))/(cal_w*(1+bet)*(1+phi*(1+1/mu_w))))*mc_w + m_w;

% 14th eq
PR(14,mc_w_p) = 1;
PR(14,w_p) = -1;
PR(14,n_p) = phi;
PR(14,lam_p) = -1;
%mc_w = w - phi*n + lam;

% 15th eq
PR(15,r_p) = 1;
PA(15,r_p) = -rho_r;
PR(15,pi_p) = -(1-rho_r)*gam_pi;
PR(15,y_p) = -(1-rho_r)*gam_y;

PRS(15,q_t) = -1;
%r = rho_r*r(-1) + (1-rho_r)*(gam_pi*pi+gam_y*y) + q;

RkP = 1/bet - (1-del);
WAP = ((alp^alp*(1-alp)^(1-alp)/(1+mu_p))*RkP^(-alp))^(1/(1-alp));
KAN = (WAP/RkP)*(alp/(1-alp));
IY = del*(KAN)^(1-alp);
RkKPY = RkP*KAN^(1-alp);

% 16th eq
PR(16,c_p) = (1/zet-IY-RkKPY);
PR(16,i_p) = (IY);
PR(16,u_p) = RkKPY;
PR(16,y_p) = -(1/zet);

PRS(16,g_t) = (1/zet);
%(1/zet-IY-RkKPY)*c + (IY)*i + RkKPY*u + (1/zet)*g = (1/zet)*y;
n_eq = 16;
%%%%%%%%%%%%%% core equations end here

%PR(17,y1_p) = 1;PA(17,y_p) = -1;
%PR(18,c1_p) = 1;PA(18,c_p) = -1;
%PR(19,i1_p) = 1;PA(19,i_p) = -1;
%PR(20,w1_p) = 1;PA(20,w_p) = -1;
%PR(21,n1_p) = 1;PA(21,n_p) = -1;
%n_eq = 21;
%%%%%%%%%%%%%% 5 lagged endogenous

G = PR(1:n_eq,1:n_eq);
H = PA(1:n_eq,1:n_eq);
FF = FU(1:n_eq,1:n_eq);

Gamma = -G;
Theta = -H;
Psi = FF;

L = FS(1:n_eq,:);
M = PRS(1:n_eq,:);

%clear PR PA FU FS PRS

s_G = size(H,1);

clear H

Xi_mat = [[Gamma  Theta] ; [eye(s_G) zeros(s_G,s_G)]];

Delta_mat = [[Psi zeros(s_G,s_G)] ; [zeros(s_G,s_G) eye(s_G)]];

m_states = s_G;
%[XI_ LA]=eig(Xi_mat,Delta_mat);

TOL = .00001;
warnings =[];
DISPLAY_IMMEDIATELY = 1;
[Xi_eigvec,Xi_eigval] = eig(Xi_mat,Delta_mat);
if rank(Xi_eigvec)<m_states,
    message = ['SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.         '
        '         Try to run your program again with DO_QZ = 1.                 '];
    if DISPLAY_IMMEDIATELY, disp(message); end;
    warnings = [warnings;message];
else
    [Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
    Xi_sortvec = Xi_eigvec(1:2*m_states,Xi_sortindex);
    Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
    Xi_select = 1 : m_states;
    if imag(Xi_sortval(m_states))~=0,
        if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) ) < TOL),
            % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE OCCURS MORE THAN ONCE!!
            % IF YOU HAVE THAT PROBLEM, PLEASE TRY MANUAL ROOT SELECTION.
            drop_index = 1;
            while (abs(imag(Xi_sortval(drop_index)))>TOL) & (drop_index < m_states),
                drop_index = drop_index + 1;
            end;
            if drop_index >= m_states,
                message = ['SOLVE.M: You are in trouble. You have complex eigenvalues, and I cannot'
                    '   find a real eigenvalue to drop to only have conjugate-complex pairs.'
                    '   Put differently: your PP matrix will contain complex numbers. Sorry!'
                    '   Try increasing the dimension of your state space. You may then get  '
                    '   sunspots, too.                                                      '];
                if DISPLAY_IMMEDIATELY, disp(message); end;
                warnings = [warnings;message];
            else
                message = ['SOLVE.M: I will drop the lowest real eigenvalue to get real PP.        '
                    '         I hope that is ok. You may have sunspots.                     '];
                if DISPLAY_IMMEDIATELY, disp(message); end;
                warnings = [warnings;message];
                Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
            end; % if drop_index >= m_states,
        end; % if (abs( Xi_sortval(m_states) - ...
    end; % if imag(Xi_sortval(m_states))~=0,
    %        if MANUAL_ROOTS,
    %          message = ['SOLVE.M: You have chosen to select roots manually.  I am crossing my   '
    %                     '         fingers that you are doing it correctly.  In particular,      '
    %                     '         you should have defined Xi_manual.  Type help solve           '
    %                     '         and inspect SOLVE.M to get further information on how to do it'];
    %          if DISPLAY_IMMEDIATELY, disp(message); end;
    %          warnings = [warnings;message];
    %          if exist('Xi_manual'),
    %             Xi_select = Xi_manual;
    %          else
    %             message = ['SOLVE.M: You have not defined Xi_manual.  Either define it or turn off '
    %                        '         the manual roots selection procedure with                     '
    %                        '         MANUAL_ROOTS = 0                                              '
    %                        '         Right now, I better let your calculations crash - sorry!      '
    %                        '         If you get results, they are based on previous calculations.  '];
    %             disp(message);
    %             warnings = [warnings;message];
    %          end; % if exist('Xi_manual'),
    %        else
    %          if max(Xi_select) < 2*m_states,
    %            if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
    %              message = ['SOLVE.M: You may be in trouble. There are stable roots NOT used for PP.'
    %                         '         I have used the smallest roots: I hope that is ok.            '
    %                         '         If not, try manually selecting your favourite roots.          '
    %                         '         For manual root selection, take a look at the file solve.m    '
    %                         '         Watch out for sunspot solutions.                              '
    %                         '         Better yet: move the time index of some endogenous variables  '
    %                         '         back by one and turn them into (predetermined) state variables'];
    %              if DISPLAY_IMMEDIATELY, disp(message); end;
    %              warnings = [warnings;message];
    %            end; % if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
    %          end; % if max(Xi_select) < 2*m_states,
    %        end; % if MANUAL_ROOTS,
    if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
        message = ['SOLVE.M: You may be in trouble.  There are unstable roots used for PP. '
            '         Keep your fingers crossed or change your model.               '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
    end; % if max(abs(Xi_sortval(Xi_select))) ...
    if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
        message = ['SOLVE.M: Your matrix PP contains a unit root. You probably do not have '
            '         a unique steady state, do you?  Should not be a problem, but  '
            '         you do not have convergence back to steady state after a shock'
            '         and you should better not trust long simulations.             '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
    end; % if abs( max(abs(Xi_sortval(Xi_select))) - 1 ...
    Lambda_mat = diag(Xi_sortval(Xi_select));
    Omega_mat  = [Xi_sortvec((m_states+1):(2*m_states),Xi_select)];
    if rank(Omega_mat)<m_states,
        message = 'SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.          ';
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
    else
        PP = Omega_mat*Lambda_mat/Omega_mat;
        PP_imag = imag(PP);
        PP = real(PP);
        if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
            message = ['SOLVE.M: PP is complex.  I proceed with the real part only.            '
                '         Hope that is ok, but you are probably really in trouble!!     '
                '         You should better check everything carefully and be           '
                '         distrustful of all results which follow now.                  '];
            if DISPLAY_IMMEDIATELY, disp(message); end;
            warnings = [warnings;message];
        end; % if sum(sum(abs(PP_imag)))
    end; % if rank(Omega_mat)<m_states,
    % End of calculating the PP matrix.  Now comes the rest.
    %calc_qrs;
end; % if rank(Xi_eigvec)<m_states,

save PP PP

% 7*16
C = [0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	0	1	0	0	0	0
    0	0	0	0	1	0	0	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	1	0	0	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	1	0	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	1	0	0	0	0	0	0
    0	0	0	0	0	0	0	0	0	0	1	0	0	0	0	0];

if restrict_v == 1;
    C = C(2:end,:);
end

n_obs = size(C,1);

H = C';
%R = zeros(size(C,1),size(C,1));
n_sh = 8;

%thetaw=theta_w;
%thetap=theta_p;
%rhop=rho_p;
%rhow=rho_w;
%rhod=rho_d;
%rhog=rho_g;
%rhoq=rho_q;

% 16*16
A = [0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
1	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	1	0	0	0	0	0	0	0	0	0	0	0	1
0	0	0	0	rho_p	-theta_p*sig_p	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	rho_w	-theta_w*sig_w	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	rho_g	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	rho_d	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	rho_q	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	0	1
0	0	0	0	0	0	0	0	0	0	0	0	1	0	0	0
0	1	0	0	0	0	0	0	0	0	0	0	0	0	0	0
0	0	0	0	0	0	0	0	0	0	0	0	0	0	1	0];
save A A
N = A;
F = A;

% (16*8)
Qmat = [1	0	0           0       0	0	0	0
        0	0	0           0       0	0	0	0
        0	1	0           0       0	0	0	0
        0	0	0           0       0	0	0	1
        0	0	1           0       0	0	0	0
        0	0	1/sig_p     0       0	0	0	0
        0	0	0           1       0	0	0	0
        0	0	0           1/sig_w	0	0	0	0
        0	0	0           0       1	0	0	0
        0	0	0           0       0	1	0	0
        0	0	0           0       0	0	1	0
        1	1	0           0       0	0	0	0
        0	0	0           0       0	0	0	1
        0	0	0           0       0	0	0	0
        0	0	0           0       0	0	0	0
        0	0	0           0       0	0	0	0];


B = Qmat;
save B B

VV = diag([sig_eps sig_v sig_p sig_w sig_g sig_d sig_q sig_ns]).^2;

D = zeros(n_obs,n_sh);
R = D*VV*D';

Q = Qmat*VV*Qmat';

%Pv = inv(eye(11^2) - kron(F,F))*reshape(Q,11^2,1)
%P = reshape(Pv,11,11);
P = .1*eye(size(F,2));
fl = 1;
t = 0;
while t < 5000
    t = t+1;
    Pnew = F*(P - P*H*inv(H'*P*H + R)*H'*P)*F' + Q;    
    fl = max(max(abs(P-Pnew)));
    P = Pnew;
    if fl <.00001;
        %t
        break
    end
end
%P
%K = F*P*H*inv(H'*P*H+R);
K_st = P*H*inv(H'*P*H+R);

%%%%%%%%%%%%%%%%%%%%%
n_exog = size(A,1);

T = (eye(n_exog)-K_st*H')*F + K_st*H'*A;% replace the dynamics for the states with kalman filter for states at time t
%%%%%%%%%%%%%%%%%%%%%

V = kron(T',FF) + kron(eye(n_exog),FF*PP+G); %
Qvec = -inv(V)*reshape(L*N+M,n_exog*m_states,1);

%check
if max(V*Qvec - (-reshape([L*N+M],n_exog*m_states,1)))>.0001;disp('Problems!!');end

RR = reshape(Qvec,m_states,n_exog);% matrix in front of the expectations of the states
save RR RR

%%%%%%%%%%%%%%%%%%%%%%%
% weight on past xi_t-1|t-1: (I-K_st*H')*F
IK_st = (eye(n_exog)-K_st*H')*F;
%IK_st
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% point C
% dynamics of the states as follow: O*[y_t x_t|t x_t s_t]' = TR*[y_t-1 x_t-1|t-1 x_t-1 s_t-1]' + SS*[shock]'

DD = n_eq+n_states+n_states+n_obs;
O = zeros(DD,DD);
O = [[eye(n_eq,n_eq) ; zeros(n_states,n_eq) ; zeros(n_states,n_eq) ; zeros(n_obs,n_eq)] [-RR ; eye(n_states) ; zeros(n_states,n_states) ; zeros(n_obs,n_states)] ...
    [zeros(n_eq,n_states) ; zeros(n_states,n_states) ; eye(n_states) ; -C] [zeros(n_eq,n_obs) ; -K_st ; zeros(n_states,n_obs) ; eye(n_obs)]];

TR = zeros(DD,DD);
TR(1:n_eq,1:n_eq) = PP;
TR(n_eq+1:n_eq+n_states,n_eq+1:n_eq+n_states) = IK_st;
TR(n_eq+n_states+1:n_eq+2*n_states,n_eq+n_states+1:n_eq+2*n_states) = A;

%%VV = diag(ones(1,7)).^2;

VAR1 = inv(O)*TR;% 51 elements
%VAR1
N = n_eq+2*n_states+n_obs;
VAR1 = [VAR1 zeros(size(VAR1,1),5)];VAR1 = [VAR1 ; zeros(5,size(VAR1,2))];% add 5 columns to have the t-1 elements
VAR1(N+1,c_p)=1;% c_t-1
VAR1(N+2,i_p)=1;% i_t-1
VAR1(N+3,y_p)=1;% y_t-1
VAR1(N+4,n_p)=1;% n_t-1
VAR1(N+5,w_p)=1;% w_t-1
%save('check2.mat','VAR1');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%load data
% variables in data: dca diad dn r pi dwa dya;
%load data pi
%data = center([dya dca diad dwa dn pi r]);
%clear pi

%n_s = size(data,2);
n_s = 7;
ZZ = zeros(n_s,size(VAR1,1));

%K = mean(data)';

ZZ(1,y_p) = 1;%Y_t
ZZ(1,N+3) = -1;%Y_t-1
ZZ(1,n_eq + dpty_t) = 1;% dpty_t

ZZ(2,c_p) = 1;%C_t
ZZ(2,N+1) = -1;%C_t-1
ZZ(2,n_eq + dpty_t) = 1;% dpty_t

ZZ(3,i_p) = 1;%i_t
ZZ(3,N+2) = -1;%i_t-1
ZZ(3,n_eq + dpty_t) = 1;% dpty_t

ZZ(4,w_p) = 1;%w_t
ZZ(4,N+5) = -1;%w_t-1
ZZ(4,n_eq + dpty_t) = 1;% dpty_t

ZZ(5,n_p) = 1;%n_t
ZZ(5,N+4) = -1;%n_t-1

ZZ(6,pi_p) = 1;%pi_t

ZZ(7,r_p) = 1;%r_t

%ZZ = ZZ*100;

DD = zeros(n_s,n_s);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

SS = [zeros(n_eq,n_sh) ; zeros(n_states,n_sh) ; Qmat*(VV^.5) ; D*(VV^.5)];

VAR_SH = inv(O)*SS;
VAR_SH = [VAR_SH ; zeros(5,n_sh)];% add 5 columns

IRF = VAR_SH;
ZIRF = ZZ*IRF;
for t = 1:50;
    IRF(:,:,t+1) = VAR1*IRF(:,:,t);
    ZIRF(:,:,t+1) = ZZ*IRF(:,:,t+1);
end

clear LIRF
for t=1:7
    if t==1|t==2|t==3|t==4|t==5;
        LIRF(t,:,:) = cumsum(ZIRF(t,:,:),3);
    else
        LIRF(t,:,:) = ZIRF(t,:,:);
    end
end


%[lnlv1,x_end,Sigma_end] = dsge_kalman(data,VAR1,VAR_SH*VAR_SH',ZZ,DD,K,zeros(size(VAR1,1),1),10*eye(size(VAR_SH*VAR_SH')));
%disp('standard, with 7 shocks')
%disp(num2str(sum(lnlv1)))
%disp(num2str((lnlv1)))

clear SER Z_SER
VV_sh = [sig_eps sig_v sig_p sig_w sig_g sig_d sig_q sig_ns]'.*randn(n_sh,1);
SS_sh = [zeros(n_eq,1) ; zeros(n_states,1) ; Qmat*(VV_sh) ; D*(VV_sh)];
clear VAR_SH 
VAR_SH = inv(O)*SS_sh;
VAR_SH = [VAR_SH ; zeros(5,1)];% add 5 columns
SHOCKS(t,:) = VV_sh';  
SER = VAR_SH;
Z_SER = ZZ*SER;

for t = 1:obs+99;
    VV_sh = [sig_eps sig_v sig_p sig_w sig_g sig_d sig_q sig_ns]'.*randn(n_sh,1);
    SS_sh = [zeros(n_eq,1) ; zeros(n_states,1) ; Qmat*(VV_sh) ; D*(VV_sh)];

    VAR_SH = inv(O)*SS_sh;
    VAR_SH = [VAR_SH ; zeros(5,1)];% add 5 columns

    SER(:,t+1) = VAR1*SER(:,t) + VAR_SH;
    Z_SER(:,t+1) = ZZ*SER(:,t+1);
    
    SHOCKS(t+1,:) = VV_sh';  
end
%size(SHOCKS)
sig_s = sqrt(sig_eps^2 + sig_v^2);

d = .2;
z2_eps = conv([0 0 0 0 1],[1 d])*(sig_v^2)/(sig_s^2);% b(L)*d(L)*(sig_v/sig_s)^2
f = .8;
z1_eps = [1 f 0 0 0 0];% f(L)
z_eps = z1_eps + z2_eps; 
Z1 = filter(z_eps,1,SHOCKS(:,1));% part of confidence driven by news shocks

z2_v = -conv([0 0 0 0 1],[1 d])*(sig_eps^2)/(sig_s^2);% b(L)*d(L)*(sig_eps/sig_s)^2
z1_v = [1 f 0 0 0 0]; % f(L)
z_v = z1_v + z2_v; 
Z2 = filter(z_eps,1,SHOCKS(:,2));% part of confidence driven by noise shocks

Conf1 = Z1 + Z2;
%size(Conf1)
Conf2 = SHOCKS(:,3:end-1)*.5*[1 1 1 1 1]';%[. .5 .6 .7 .8]';% part of confidence driven by the "known" shocks (this is zero on impact)
%size(Conf2)
%pause

%other = SHOCKS(:,3:end)*[.08 .04 .04 .06 .08]';other = other(100:end-1);

SER = SER(:,101:end)';
Z_SER = Z_SER(:,101:end)';
SHOCKS = SHOCKS(101:end,:);

Conf = Conf1(101:end,:) + Conf2(100:end-1,:);
save Conf Conf

clear CUM_SER
for t=1:7
    if t==1|t==2|t==3|t==4|t==5;
        CUM_SER(:,t) = cumsum(Z_SER(:,t));
    else
        CUM_SER(:,t) = Z_SER(:,t);
    end
end

%CUM_SER = [SER(:,20)+other CUM_SER];% first element of CUM_SER is the technology process
CUM_SER = [SER(:,20) CUM_SER];% first element of CUM_SER is the technology process

%lab{1}='Y';lab{2}='C';lab{3}='I';lab{4}='W';lab{5}='N';lab{6}='\pi';lab{7}='r';lab{8}='tech';
%for t=1:8
%    subplot(2,4,t);plot(CUM_SER(:,t));title(lab{t});
%end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



