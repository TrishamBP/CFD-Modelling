% One-dimensional unsteady diffusion by the FTCS and Crank-Nicolson schemes

clear all;
close all;
clc;

% scheme=1 <===> FTCS Explicit
% scheme=2 <===> FTCS Implicit
% scheme=3 <===> Crank-Nicolson
scheme=3;

time=0.0;
timefinal=0.18;%1.04;
nu=0.000217;
Uo=40;

dt=0.002;
H=0.04;
dy=0.001;
Y=0:dy:H;
NY=length(Y);
BB=(nu*dt)/(dy^2);

% Numerical solution at n step
fo(1:NY)=0;
% Numerical solution at n+1
f=fo;
% Exact solution
fex(1:NY)=0;

figure(1);
axis([0 40 0 H])

if(scheme==1)
    disp('FTCS Explicit')
    %%% FTCS Explicit
    
    tic
    
    while(time<=timefinal)
        time=time+dt;
        plot(f,Y,'b',fex,Y,'or')
        title(time)
        pause(0.01)
        
        % Numerical solution
        f(1)=Uo;
        f(NY)=0;
        fo=f;
        
        for i=2:NY-1
            f(i)=fo(i)+BB*(fo(i+1)-2*fo(i)+fo(i-1));
        end
        
        % Exact Solution
        eta=Y./(2*sqrt(nu*time));
        eta1=H/(2*sqrt(nu*time));
        
        erfcsum=erfc(eta);
        
        for k=1:100
            erfcsum=erfcsum-erfc(2*k*eta1-eta)+erfc(2*k*eta1+eta);
        end
        
        fex=Uo.*erfcsum;
        
    end;
    toc
end

if (scheme==2)
    disp('FTCS Implicit')
    
    rhs(2:NY-1)=0;
    %%% FTCS Implicit
    
    tic
    
    %build the matrix
    AA(2:NY-1,2:NY-1)=0;
    AA(2,2)=2*BB+1;
    AA(2,3)=-BB;
    AA(NY-1,NY-1)=2*BB+1;
    AA(NY-1,NY-2)=-BB;
    
    for jj=3:NY-2
        AA(jj,jj)=2*BB+1;
        AA(jj,jj+1)=-BB;
        AA(jj,jj-1)=-BB;
    end
    
    while(time<=timefinal)
        time=time+dt;
        plot(fex,Y,'b',f,Y,'or')
        title(time)
        pause(0.01)
        
        % Numerical solution
        f(1)=Uo;
        f(NY)=0;
        fo=f;
        
        %build the matrix and rhs
        rhs(3:NY-2)=fo(3:NY-2);
        rhs(2)=fo(2)+BB*fo(1);
        rhs(NY-1)=fo(NY-1)+BB*fo(NY);
        
        % Solve system of linear equations using matlab intrincics
        
        f(2:NY-1)=AA(2:NY-1,2:NY-1)\rhs(2:NY-1)';
        
        % Exact Solution
        eta=Y./(2*sqrt(nu*time));
        eta1=H/(2*sqrt(nu*time));
        
        erfcsum=erfc(eta);
        
        for k=1:100
            erfcsum=erfcsum-erfc(2*k*eta1-eta)+erfc(2*k*eta1+eta);
        end
        
        fex=Uo.*erfcsum;
        
    end;
    toc
end


if (scheme==3)
    disp('Crank-Nicolson')
    
    rhs(2:NY-1)=0;
    %%% Crank-Nicolson
    
    tic
    
    %build the matrix
    AA(2:NY-1,2:NY-1)=0;
    AA(2,2)=BB+1;
    AA(2,3)=-0.5*BB;
    AA(NY-1,NY-1)=BB+1;
    AA(NY-1,NY-2)=-0.5*BB;
    
    for jj=3:NY-2
        AA(jj,jj)=BB+1;
        AA(jj,jj+1)=-0.5*BB;
        AA(jj,jj-1)=-0.5*BB;
    end
    
    while(time<=timefinal)
        time=time+dt;
        plot(fex,Y,'b',f,Y,'or')
        title(time)
        pause(0.01)
        
        % Numerical solution
        f(1)=Uo;
        f(NY)=0;
        fo=f;
        
        % build the rhs
        for i=3:NY-2
            rhs(i)=fo(i)+0.5*BB*(fo(i+1)-2*fo(i)+fo(i-1));
        end
        
        rhs(2)=BB*fo(1)+(1.0-BB)*fo(2)+0.5*BB*fo(3);
        rhs(NY-1)=BB*fo(NY)+(1.0-BB)*fo(NY-1)+0.5*BB*fo(NY-2);
        
        % Solve system of linear equations using matlab intrincics
        f(2:NY-1)=AA(2:NY-1,2:NY-1)\rhs(2:NY-1)';
        
        % Exact Solution
        eta=Y./(2*sqrt(nu*time));
        eta1=H/(2*sqrt(nu*time));
        
        erfcsum=erfc(eta);
        
        for k=1:100
            erfcsum=erfcsum-erfc(2*k*eta1-eta)+erfc(2*k*eta1+eta);
        end
        
        fex=Uo.*erfcsum;
        
    end;
    toc
end
