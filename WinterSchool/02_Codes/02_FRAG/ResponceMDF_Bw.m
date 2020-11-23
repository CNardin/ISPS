function [HistVar] = ResponceMDF_Bw(Mat,Var)
%% General Parameters
nelem=Mat.nelem;
Af=Mat.Af;
%% Force parameters 
Pext=Mat.Fe;
%% Time variables 
dt=Mat.dt;
t=Mat.t;
%% Define Parameters of Newmark scheme 
gamma=Mat.gammaN;
beta=Mat.betaN;
a1 = 1.0/(beta*dt*dt);
a2 = -1.0/(beta*dt*dt);
a3 = -1.0/(beta*dt);
a4 = 1.0 - 1.0/(2.0*beta);
a5 = gamma/(beta*dt);
a6 = -gamma/(beta*dt);
a7 = 1.0 - gamma/beta;
a8 = dt*(1.0 - gamma/(2.0*beta));

%%
HistVar.e=zeros(nelem,numel(t));
HistVar.z=zeros(nelem,numel(t));
%HistVar.k=zeros(nelem,numel(t));
HistVar.u=zeros(nelem,numel(t));
HistVar.Q=zeros(nelem,numel(t));
HistVar.du=zeros(nelem,numel(t));
HistVar.ddu=zeros(nelem,numel(t));
HistVar.eps=zeros(nelem,numel(t));
%HistVar.k(1)=Mat.ko;
HistVar.grad=zeros(nelem,numel(t));


MatCur.Z(1:nelem,1)=0;
MatCur.E(1:nelem,1)=0;
MatCur.Q(1:nelem,1)=0;


%% Initial Condition 

u=Mat.u;
du=Mat.du;
ddu=Mat.ddu;
%% Initialization 
Var.Z(1:nelem,1)=0;
Var.E(1:nelem,1)=0;

%% Loop

for i=1:numel(t)
 %   disp(i);
    Res(1:nelem,1)=1; %spurius value
    d(1:nelem,1)=1; %spurius value 
    uold=u;
    duold=du;
    dduold=ddu;
    Ptilda0=Pext(:,i)-Mat.M*(a2*uold + a3*duold + a4*dduold)...
           -Mat.C*( a6*uold + a7*duold + a8*dduold);
    Res0=((Mat.M*a1+Mat.C*a5+Mat.ke)*u-Ptilda0);
    Var.DEPSI=Af*(u-uold);
    Var.EPS=Af*u;  
    for j=1:nelem; 
        Mat.alpha=Mat.ALPHA(j);
        Mat.ko=Mat.KO(j);
        Mat.n=Mat.N(j);
        Mat.Ao=Mat.AO(j);
        Mat.deltaA=Mat.DELTAA(j);
        Mat.deltav=Mat.DELTAv(j);  
        Mat.deltaeta=Mat.DELTAeta(j);
        %
        Var.z=Var.Z(j);
        Var.e=Var.E(j);
        Var.depsi=Var.DEPSI(j);
        Var.eps=Var.EPS(j);
        [MatCurs] = BoucwenBN(Mat,Var);
        KsTemp(j,j)=MatCurs.k;
        MatCur.Q(j,1)=MatCurs.sig;
        MatCur.Z(j,1)=MatCurs.z;
        MatCur.E(j,1)=MatCurs.e;
    end
    CurK=Af'*KsTemp*Af;
    Ktilda=(Mat.M*a1+Mat.C*a5+CurK);
    d1=Ktilda\(Ptilda0-(Mat.M*a1+Mat.C*a5)*u-Af'*MatCur.Q);
    u=u+d1;
           Var.DEPSI=Af*(u-uold);
        Var.EPS=Af*u;
        for j=1:nelem; 
            Mat.alpha=Mat.ALPHA(j);
            Mat.ko=Mat.KO(j);
            Mat.n=Mat.N(j);
            Mat.Ao=Mat.AO(j);
            Mat.deltaA=Mat.DELTAA(j);
            Mat.deltav=Mat.DELTAv(j);  
            Mat.deltaeta=Mat.DELTAeta(j);
            %
            Var.z=Var.Z(j);
            Var.e=Var.E(j);
            Var.depsi=Var.DEPSI(j);
            Var.eps=Var.EPS(j);
            [MatCurs] = BoucwenBN(Mat,Var);
            KsTemp(j,j)=MatCurs.k;
            MatCur.Q(j,1)=MatCurs.sig;
            MatCur.Z(j,1)=MatCurs.z;
            MatCur.E(j,1)=MatCurs.e;
        end    
%    
    if abs(d1'*Res0)~=0;
    while (abs(d'*Res)/abs(d1'*Res0))>10^-10
        Ptilda=Pext(:,i)-Mat.M*(a2*uold + a3*duold + a4*dduold)...
           -Mat.C*( a6*uold + a7*duold + a8*dduold);
        Ktilda=(Mat.M*a1+Mat.C*a5+CurK);
        d=Ktilda^-1*(Ptilda-(Mat.M*a1+Mat.C*a5)*u-Af'*MatCur.Q);
        u=u+d;
        Var.DEPSI=Af*(u-uold);
        Var.EPS=Af*u;

        for j=1:nelem;
            % material 
            Mat.alpha=Mat.ALPHA(j);
            Mat.ko=Mat.KO(j);
            Mat.n=Mat.N(j);
            Mat.Ao=Mat.AO(j);
            Mat.deltaA=Mat.DELTAA(j);
            Mat.deltav=Mat.DELTAv(j);  
            Mat.deltaeta=Mat.DELTAeta(j);
            %
            Var.z=Var.Z(j);
            Var.e=Var.E(j);
            Var.depsi=Var.DEPSI(j);
            Var.eps=Var.EPS(j);
            [MatCurs] = BoucwenBN(Mat,Var);
            KsTemp(j,j)=MatCurs.k;
            MatCur.Q(j,1)=MatCurs.sig;
            MatCur.Z(j,1)=MatCurs.z;
            MatCur.E(j,1)=MatCurs.e;
        end
        CurK=Af'*KsTemp*Af;
        Res=(Ptilda-(Mat.M*a1+Mat.C*a5)*u-Af'*MatCur.Q);
    end
    end
    du=a5*u+a6*uold+a7*duold+a8*dduold;
    ddu=a1*u+a2*uold+a3*duold+a4*dduold;
    %% Update the values 
    Var.EPS=Af*u;
    Var.Z=MatCur.Z;
    Var.e=MatCur.E;
    Var.Q=MatCur.Q;
    Var.k=CurK;
    Var.dEPSI=Af*(u-uold);


%% Store Values   
if i==1
    HistVar.u(:,i)=u;
    HistVar.z(:,i)=Var.Z;
    HistVar.e(:,i)=MatCur.E;
    HistVar.Q(:,i)=MatCur.Q;
    HistVar.k{i}=CurK;
    HistVar.du(:,i)=du;
    HistVar.ddu(:,i)=ddu;
    HistVar.eps(:,i)=Af*u;
else
    HistVar.z(:,i)=Var.Z;
    HistVar.e(:,i)=MatCur.E;
    HistVar.Q(:,i)=MatCur.Q;
    HistVar.k{i}=CurK;
    HistVar.u(:,i)=u;
    HistVar.du(:,i)=du;
    HistVar.ddu(:,i)=ddu;
    HistVar.eps(:,i)=Af*u;
end 
end
end