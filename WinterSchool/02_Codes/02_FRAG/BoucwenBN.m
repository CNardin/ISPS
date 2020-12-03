%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%________________U N I V E R S I T Y    O F   C A L I F O R N I A _____________
%_____________________________B E R K E L E Y__________________________________
%_________________________________CE 234_______________________________________
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function written by           marco broccardo                  April  22th 2011 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%===============================================================================
%               marco broccardo (marcobroccardo@berkeley.edu)
%  Department of Civil and Environmental Engineering, UC Berkeley
%  Copyright(c) The Regents of the University of California. All Rights Reserved.
%===============================================================================
%===============================================================================
%====================== Uniaxial Bouc-Wen Material==========================

%% Material Parameters
function [MatCur] = BoucwenBN(Mat,Var) 
alpha=Mat.alpha;
ko=Mat.ko;
n=Mat.n;
gamma=Mat.gamma;
beta=Mat.beta;
Ao=Mat.Ao;
deltaA=Mat.deltaA;
deltav=Mat.deltav;
deltaeta=Mat.deltaeta;
%%
tol=10^-12;
znew=1;
zold=0;
z=Var.z;
    while abs(znew-zold)>tol;
          abs(znew-zold);
          %evaluate the function 
          e=Var.e+(1-alpha)*ko*(Var.depsi)*z;
          A=Ao-deltaA*e;
          ni=1+deltav*e;
          eta=1+deltaeta*e;
          PSI=gamma+beta*sign((Var.depsi)*z);
          PHI=A-abs(z)^n*PSI*ni;
          f=z-Var.z-PHI/eta*(Var.depsi);
          %evaluate the function derivatives 
          ep=(1-alpha)*ko*(Var.depsi);
          Ap=-deltaA*ep;
          nip=deltav*ep;
          etap=deltaeta*ep;
          PHIp=Ap-n*abs(z)^(n-1)*sign(z)*PSI*ni-abs(z)^n*PSI*nip;
          fp=1-(PHIp*eta-PHI*etap)/eta^2*(Var.depsi);
          % New trial value in the Newton scheme 
          znew=z-f/fp;
          % Update 
          zold=z;
          z=znew;   
    end   
    MatCur.sig=alpha*ko*Var.eps+(1-alpha)*ko*z;
    %Consistent tangent
    e=Var.e+(1-alpha)*ko*(Var.depsi)*z;
    A=Ao-deltaA*e;
    ni=1+deltav*e;
    eta=1+deltaeta*e;
    PSI=gamma+beta*sign((Var.depsi)*z);
    PHI=A-(abs(z))^n*PSI*ni;    
    b1=(1-alpha)*ko*z;
    b2=(1-alpha)*ko*((Var.depsi));
    b3=((Var.depsi))/eta;
    b4=-b3*deltaA*b1-b3*abs(z)^n*PSI*deltav*b1-PHI/(eta^2)*(Var.depsi)*deltaeta*b1...
        +PHI/eta;
    b5=1+b3*deltaA*b2+b3*n*abs(z)^(n-1)*sign(z)*PSI*ni+b3*abs(z)^n*PSI*deltav*b2...
        +PHI/(eta^2)*(Var.depsi)*deltaeta*b2;
    
    MatCur.k=alpha*ko+(1-alpha)*ko*b4/b5;
    % Stack current information
    MatCur.z=z;
    MatCur.e=e;
end

