function y= penalty_mp_shock(q1)
global IR
global ssigma;
global nper;


index=[4 5];
x=zeros(nper,size(index,2));

for i=1:size(index,2);
    for h=1:nper
        x(h,i)=-(1/ssigma(index(i)))*squeeze(IR(h,index(i),:))'*q1;
		if i ==2
			x(h,i) = -x(h,i);
		end
    end
end

y=0;
for i=1:size(index,2);
    for h=1:nper

    if (x(h,i)>0)
        y=100*x(h,i)+ y;
    else
        y=x(h,i) + y;
    end
    
    end
end



end




