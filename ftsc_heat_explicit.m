function h_new=ftsc_heat_explicit(x_num,x,t,dt,cfl,rhs,bc,h)

h_new=zeros(x_num,1);

L=1:x_num-2;
C=2:x_num-1;
R=3:x_num;

f=rhs(x_num,x,t);

h_new(C)=h(C)+cfl*(h(L)-2.0*h(C)+h(R))+dt*f(C);

h_new=bc(x_num,x,t+dt,h_new);
