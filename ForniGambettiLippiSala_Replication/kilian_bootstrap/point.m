% POINT.M 
% Lutz Kilian
% University of Michigan
% April 1997


function [CI]=point(y)

global h p q

p=4;		% VAR lag order
h=40;		% Impulse response horizon

[t,q]=size(y);
y=detrend(y,0);
[A,SIGMA,U,V]=olsvarc(y,p);						% VAR with intercept
if ~ any(abs(eig(A))>=1)
	[A]=asybc(A,SIGMA,t,p);
end;
[CI]=boot(A,U,y,V);
