function G=calculateG(Fault,XY,los)
% calculate Green function 
Xe=XY(:,1);
Yn=XY(:,2);
xs=Fault(:,1);
ys=Fault(:,2);
zs=Fault(:,3);
strike=Fault(:,16);
dip=Fault(:,17);
len=Fault(:,18);
width=Fault(:,19);



for i = 1:length(xs)
    for j = 1:length(Xe)
        xrs(i,j) =Xe(j)-ys(i);%x distance, in km
        yrs(i,j) =Yn(j)-xs(i);%y distance, in km
        zrs(i,j) = -zs(i);%z distance, in km
    end
end

[G] = okada_InSAR(xrs,yrs,zrs,strike,dip,len,width,los(:,1),los(:,2),los(:,3));

