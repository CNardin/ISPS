

%% Geometry 
Mat.nelem=5;        % number of elements 
Mat.NDOF=5;         % number of degree of freedom
Mat.Af= 5 ;         % compatibility matrix  

%% Mass matrix
me=1*10^6;          % elementary mass 
Mat.M=me;           % mass matrix
Mat.M(1,1)=5.04*10^6; %[kg]
Mat.M(2,2)=3.22*10^6; 
Mat.M(3,3)=3.22*10^6;
Mat.M(4,4)=3.02*10^6;
Mat.M(5,5)=1*2.55*10^6; 

%% Stiffness matrix  
ke=2.8*10^8;        % elementary stiffness 
Mat.ke=ke;          % stiffness matrix 
ks(1,1)=1*4.3790*10^7;%[kN/m]
ks(2,2)=1*4.0675*10^7;
ks(3,3)=0.9*3.9560*10^7;
ks(4,4)=0.7*3.9150*10^7;
ks(5,5)=0.5*3.9575*10^7;
Mat.ke=Mat.Af'*ks*Mat.Af;

%% Damping matrix 
we=sqrt(ke/me);     % characteristic frequency 
ce=2*we*0.05*me;    % damping of the element 
cs(1,1) = 2*0.05*sqrt(ks(1,1)*Mat.M(1,1));
cs(2,2) = 2*0.05*sqrt(ks(2,2)*Mat.M(2,2));
cs(3,3) = 2*0.05*sqrt(ks(3,3)*Mat.M(3,3));
cs(4,4) = 2*0.05*sqrt(ks(4,4)*Mat.M(4,4));
cs(5,5) = 2*0.05*sqrt(ks(5,5)*Mat.M(5,5));
Mat.C=Mat.Af'*cs*Mat.Af; % damping matrix 
       
%% Check 
[V,D] = eig(Mat.ke,Mat.M);
%V
TT=2*pi./sqrt(D);   % structural period
FF=1./TT;           % structural frquency 

%% Material Parameters  
alphaaa=.01;         % Post yielding stiffness 
%threshold=0.12;
uy=0.01;            % yielding point

%% Define Parameters of Newmark scheme 
Mat.gammaN=1/2;
Mat.betaN=1/6;


Mat.u(1:Mat.NDOF,1)=0;
Mat.du(1:Mat.NDOF,1)=0;
Mat.ddu(1:Mat.NDOF,1)=0;
%
Mat.r(1:Mat.NDOF)=1;
% Mat.r = Mat.r';

TT
