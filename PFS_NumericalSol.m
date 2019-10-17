% Numerical Solution of the Porous-Fisher-Stefan model.
% Inputs: Stefan constant (kappa), final time (tf)
% Output: Speed of travelling wave (Ck)

function Ck=PFS_NumericalSol(kappa,tf)
global a b c d N

dt=1e-2;

tol=1e-7;
picmax=1000;
N=5e4+1;
dxi=1/(N-1);

%initial conditions
L0=10.0;
L=L0;
LUpdate=L;
U=zeros(N,1)';
UUpdate=U;

a=zeros(N,1)';
b=zeros(N,1)';
c=zeros(N,1)';
d=zeros(N,1)';

xi=0:dxi:1;

% Initial condition of U (u0) is approximation of travelling wave solution,
% valid when kappa is large

Z=(kappa+3)*sqrt(2)/(kappa+4)*( log(1-xi) + log(1+(kappa+3)*xi)/(kappa+3));
Z(end)=-1e5;
U=fliplr(interp1(-Z/L0,xi,xi,'makima')).^2;
u0=sqrt(U);


t=0.0;
MF=L; %MF records L(t)
time = 0:dt:tf;
maxsteps=length(time)-1;

for i=1:maxsteps
    dt=time(i+1)-time(i);
    t=t+dt;
    PU=U;
    PL=L;
    
    
    for pic=1:picmax
        L=LUpdate;
        U=UUpdate;
        
        
        a(1)=0.0;
        b(1)=-1.0;
        c(1)=1.0;
        d(1)=0.0;
        
        a(N)=0;
        b(N)=1.0;
        c(N)=0;
        d(N)=0.0;
        
        
        a(2:N-1)=sqrt(abs(U(2:N-1))) -...
            (L*dxi)*xi(2:N-1)*(L-PL)/(2.0*dt);
        
        b(2:N-1)=-2 .*sqrt(abs(U(2:N-1))) -...
            (L*dxi).^2./dt+ 2*(L*dxi).^2*(1.0-sqrt(abs(U(2:N-1))));
        
        c(2:N-1)=sqrt(abs(U(2:N-1))) +...
            (L*dxi).*xi(2:N-1)*(L-PL)/(2.0*dt);
        
        d(2:N-1)=-1.*PU(2:N-1)/dt*(L*dxi).^2;
        
        %Thomas Algorithm to update U
        UUpdate = thomas(N,a,b,c,d);
        
        %Update the moving boundary L(t)
        LUpdate = sqrt(PL^2-kappa*dt*...
            (3*UUpdate(N)-4*UUpdate(N-1)+UUpdate(N-2))/(2*dxi));
        
        
        if norm(U-UUpdate,Inf) < tol && abs(L-LUpdate)<tol
            if mod(i,100)==0
                fprintf('Time is %d\n',t);
                fprintf('Iterations %d\n',pic);
            end
            break
        end
        if pic==picmax
            error('ERROR: Max Picard Iteration Reached')
            
        end
    end
    
    PU=UUpdate;
    PL=LUpdate;
    MF(i+1)=LUpdate;
end

x=LUpdate*xi;

figure(900)
hold off
plot(x-L0,sqrt(abs(UUpdate)),'k','LineWidth',2)
hold on

plot([0 linspace(x(end)-2*L0,x(end)-L0,N)],[1 u0],'r--','LineWidth',1.5);

axis([0 30 0 1])
set(gca,'FontSize',18)
ylabel('u(x-ct)','FontSize',18)
xlabel('x','FontSize',18)


M=length(MF);
disp(time(M))


figure(901)
%solution of the moving boundary L(t)
plot(time(1:M),MF)

%speed of the travelling wave c = L'(t). This should be roughly constant.
plot(time(2:M),diff(MF)/dt); grid on;

%slope of best-fit line near tf determines travelling wave speed
A=polyfit(time(M-10:M),MF(M-10:M),1);
Ck=A(1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Subroutines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Thomas Algorithm
function x = thomas(N,a,b,c,d)
bb=b;
dd=d;

for i=2:N
    ff=a(i)/bb(i-1);
    bb(i)=bb(i)-c(i-1)*ff;
    dd(i)=dd(i)-dd(i-1)*ff;
end
x(N)=dd(N)/bb(N);

for i=1:N-1
    j=N-i;
    x(j)=(dd(j)-c(j)*x(j+1))/bb(j);
end
end


