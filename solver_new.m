clc;

m = 0;
time = tmax;
length = L; 

dx1 = dx ; % [cm]
dt1 = dt ; % [min]

x = (1:dx1:length);
t = (0:dt1:time);
sol = pdepe(m,@pdepedx1_new,@pdepeic_new,@pdepebc_new,x,t);

plot(t,sol)
ylim([7999, 8000]);
