

%% Geometry 
Mat.nelem=1; % number of elements 
Mat.NDOF=1;  % number of degree of freedom
Mat.Af= 1 ;  % compatibility matrix (Af matrix) Filippou lecture notes 

%% Mass matrix 
me=1*10^6; %elementary mass 
Mat.M=me;  % mass matrix
%% Stiffness matrix  
ke=2.8*10^8; % elementary stiffnesses 
Mat.ke=ke;   % stiffness matrix 

%% Damping matrix 
we=sqrt(ke/me); % important freqeuncy 
ce=2*we*0.05*me; % damping of the element 
Mat.C=ce; % damping matrix 
%% Check 
[V,D] = eig(Mat.ke,Mat.M);
%V
TT=2*pi./sqrt(D); % dstructural period
FF=1./TT; % strutural frquency 

%% Material Parameters  
alphaa=.01;   % Post yelding stiffness 
%threshold=0.12;
uy=0.04; % yielding point

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Buc Wen Model parameters%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

n=10;
%
Mat.ALPHA(1)=alphaa;
%
Mat.KO(1)=ke;
%
Mat.N(1:Mat.nelem)=n;
%
Mat.gamma(1)=1/(2*(uy)^n);
%
Mat.beta(1:Mat.nelem)=Mat.gamma(1);
Mat.AO(1:Mat.nelem)=1;
Mat.DELTAA(1:Mat.nelem)=0;
Mat.DELTAv(1:Mat.nelem)=0;
Mat.DELTAeta(1:Mat.nelem)=0;



%% Define Parameters of Newmark scheme 
Mat.gammaN=1/2;
Mat.betaN=1/6;


Mat.u(1:Mat.NDOF,1)=0;
Mat.du(1:Mat.NDOF,1)=0;
Mat.ddu(1:Mat.NDOF,1)=0;
%
Mat.r=1;


