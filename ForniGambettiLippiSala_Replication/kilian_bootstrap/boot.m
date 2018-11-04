% BOOT.M 
% Lutz Kilian
% University of Michigan
% April 1997

function [CI]=boot(A,U,y,V)

global p h 

nrep=2000; % for 90% interval

p=4;
h=40;

[t,q]=size(y);				
y=y';
Y=y(:,p:t);	
for i=1:p-1
 	Y=[Y; y(:,p-i:t-i)];		
end;

Ur=zeros(q*p,t-p);   
Yr=zeros(q*p,t-p+1); 
IRFrmat=zeros(nrep,q^2*(h+1));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  start of bootstrap simulation                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Create nboot bootstrap replications of pseudo data

for j=1:nrep
   
	pos=fix(rand(1,1)*(t-p+1))+1;
	Yr(:,1)=Y(:,pos);

	index=fix(rand(1,t-p)*(t-p))+1;
	Ur(:,2:t-p+1)=U(:,index);	

	for i=2:t-p+1
		Yr(:,i)= V + A*Yr(:,i-1)+Ur(:,i); 
	end;

	yr=[Yr(1:q,:)];
	for i=2:p
		yr=[Yr((i-1)*q+1:i*q,1) yr];
   end;
   yr=yr';
   yr=detrend(yr,0);
   
	pr=p;

   [Ar,SIGMAr]=olsvarc(yr,pr);

	if ~ any(abs(eig(Ar))>=1)
		[Ar]=asybc(Ar,SIGMAr,t,pr);
	end;

	[IRFr]=irfvar(Ar,SIGMAr(1:q,1:q),pr);
    size(IRFr)
	IRFrmat(j,:)=vec(IRFr)';
    CI(:,:,j) = IRFr; 
end;   

% Calculate 90 perccent interval endpoints
%CI=prctile(IRFrmat,[5 95]);
