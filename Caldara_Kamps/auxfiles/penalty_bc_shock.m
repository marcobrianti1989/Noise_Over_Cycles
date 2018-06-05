function y= penalty_bc_shock(q1)

global IR
global ssigma;
global nper;


index=[1 3];
x=zeros(nper,size(index,2));

for i=1:size(index,2);
    for h=1:nper
        x(h,i)=-(1/ssigma(index(i)))*squeeze(IR(h,index(i),:))'*q1;
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




