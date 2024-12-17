function Fault = simulatefault(Length,Width,EX,EY,EZ,Strike,Dip,f)
dip = Dip*pi/180;
arfa = (Strike)*pi/180;
ast = f(1);
adi = f(2);
pat_l = linspace(EY,Length+EY,ast);
pat_w = linspace(EX,-Width+EX,adi);
F=[];
for j=1:adi-1
    for i=1:ast-1
        F0=[pat_l(i),pat_w(j),pat_l(i),pat_w(j+1),pat_l(i+1),pat_w(j+1),pat_l(i+1),pat_w(j) Strike Dip];
        F=[F;F0];
    end
end
lowl = F(:,1:2);
upl = F(:,3:4);
upr = F(:,5:6);
lowr = F(:,7:8);
strike = F(:,9);
dip0 = F(:,10);
%Get along strike and dip lengths
strike_length = Length/(ast-1)*ones(size(strike));
dip_length = Width/(adi-1)*ones(size(dip0));
area = (strike_length.*dip_length);
R1 = [cos(arfa) sin(arfa)*cos(dip);sin(arfa) -cos(arfa)*cos(dip);0 sin(dip)];
lowlXYZ = [R1*[lowl(:,1)'-EY;lowl(:,2)'-EX]]';
uplXYZ = [R1*[upl(:,1)'-EY;upl(:,2)'-EX]]';
uprXYZ = [R1*[upr(:,1)'-EY;upr(:,2)'-EX]]';
lowrXYZ = [R1*[lowr(:,1)'-EY;lowr(:,2)'-EX]]';
xf1=lowlXYZ(:,1)+EY;xf2=uplXYZ(:,1)+EY;xf3=uprXYZ(:,1)+EY;xf4=lowrXYZ (:,1)+EY;
yf1=lowlXYZ(:,2)+EX;yf2=uplXYZ(:,2)+EX;yf3=uprXYZ(:,2)+EX;yf4=lowrXYZ (:,2)+EX;
zf1=lowlXYZ(:,3)+EZ;zf2=uplXYZ(:,3)+EZ;zf3=uprXYZ(:,3)+EZ;zf4=lowrXYZ (:,3)+EZ;
xc=(xf1+xf3)/2;yc=(yf1+yf3)/2;mzd=(zf1+zf3)/2;
Fault = [xc yc mzd xf1 xf2 xf3 xf4 yf1 yf2 yf3 yf4 zf1 zf2 zf3 zf4 strike dip0 strike_length dip_length area];