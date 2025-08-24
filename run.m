

clear,clc,close all


global Fw Fe Fs Fn DF aW aE aS aN aP bP dU dV


% Geometry
H = 0.01;       
L = 10*H;       


Nx = 200;       
dx = L/Nx;     
Ny = 40;        
dy = H/Ny;     

dz = 0.01;      

x  = dx/2:dx:L-dx/2;    
xu = 0:dx:L;            
y  = dy/2:dy:H-dy/2;    
yv = 0:dy:H;           

iu = 1:Nx+1; Ju = 2:Ny+1;   
Iv = 2:Nx+1; jv = 1:Ny+1;   
Ip = 2:Nx+1; Jp = 2:Ny+1;  
iF = 2:Nx;   jF = 2:Ny;     


rho   = 1.2;            
mu    = 1.8e-5;        
nu    = mu/rho;         
kt    = 0.025;          
cp    = 1006;           
alpha = kt/(rho*cp);    
Pr    = nu/alpha;      


Re = 100;              
U  = Re*nu/(2*H);       
Ti = 20;                
Tw = 100;               
qw = 100;               

BC_N = 1;              
BC_S = 1;               


alphaU  = 0.3;  
alphaP  = 0.2;  
NmaxSIM = 1e+4; 
NmaxGSI = 1e+1; 
err     = 1e-5; 
div     = 1e+1; 


u      = zeros(Nx+1,Ny+2); v      = zeros(Nx+2,Ny+1);
uStar  = zeros(Nx+1,Ny+2); vStar  = zeros(Nx+2,Ny+1);
uPrime = zeros(Nx+1,Ny+2); vPrime = zeros(Nx+2,Ny+1);
dU     = zeros(Nx+1,Ny+2); dV     = zeros(Nx+2,Ny+1);

T      = zeros(Nx+2,Ny+2);
p      = zeros(Nx+2,Ny+2); 
pPrime = zeros(Nx+2,Ny+2);

Fe = zeros(Nx+1,Ny+1); Fw = zeros(Nx+1,Ny+1);   
Fn = zeros(Nx+1,Ny+1); Fs = zeros(Nx+1,Ny+1);
DF = zeros(Nx+1,Ny+1);

aE = zeros(Nx+1,Ny+1); aW = zeros(Nx+1,Ny+1);   
aN = zeros(Nx+1,Ny+1); aS = zeros(Nx+1,Ny+1);   
aP = zeros(Nx+1,Ny+1); bP = zeros(Nx+1,Ny+1);

ures  = zeros(NmaxSIM,1);  
vres  = zeros(NmaxSIM,1);   
pres  = zeros(NmaxSIM,1);   


u(:,Ju) = U;


p1  = 12*mu*U*L/(2*H)^2;    
p(Ip,Jp) = ones(Nx,Ny).*linspace(p1,0,Nx)';

T(:,Jp) = Ti; T(:,1) = Tw; T(:,Ny+2) = Tw;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


Dx = (mu/dx)*dy*dz;
Dy = (mu/dy)*dx*dz;

for n = 1:NmaxSIM

    uOld  = u;
    vOld  = v;
    pStar = p;
    
    

    FVM_u(Nx,Ny,dx,dy,dz,rho,Dx,Dy,iF,Ju,alphaU,uOld,vOld,pStar,BC_S);
    
    [uStar,ures(n)] = FVM_GS_ext_mesh(Nx,Ny+1,alphaU,NmaxGSI,err,uOld);


    FVM_v(Nx,Ny,dx,dy,dz,rho,Dx,Dy,Iv,jF,alphaU,u,v,pStar)
    
    [vStar,vres(n)] = FVM_GS_ext_mesh(Nx+1,Ny,alphaU,NmaxGSI,err,vOld);
    

    FVM_pcorr(Nx,Ny,dx,dy,dz,rho,Ip,Jp,uStar,vStar)  

    pPrime(:,:) = 0;
    [pPrime,pres(n)] = FVM_GS_ext_mesh(Nx+1,Ny+1,1,NmaxGSI,err,pPrime);  

     

    p(Ip,Jp)      = pStar(Ip,Jp) + pPrime(Ip,Jp)*alphaP;
    
    uPrime(iF,Ju) = dU(iF,Ju).*(pPrime(iF,Ju) - pPrime(iF+1,Ju));                      
    u(iF,Ju)      = uStar(iF,Ju) + uPrime(iF,Ju);
                 
    vPrime(Iv,jF) = dV(Iv,jF).*(pPrime(Iv,jF) - pPrime(Iv,jF+1));                      
    v(Iv,jF)      = vStar(Iv,jF) + vPrime(Iv,jF);
    
    
    
    if n > 10        
        fprintf('n = %5.0f, u = %6.2e, v = %6.2e, p = %6.2e \n',...
                 n,ures(n),vres(n),pres(n))
        cTest = max([ures(n),vres(n)]);
        if cTest < err
            break; 
        elseif cTest > div || isnan(cTest)
            fprintf('Residuals are too high.')
            break;
        end
    end
    
    u(Nx+1,:) = u(Nx,:);      
    v(Nx+2,:) = v(Nx+1,:);
    
end

 

% Setup
FVM_phi(Nx,Ny,dx,dy,dz,rho,kt/cp,qw/cp,Ip,Jp,u,v,BC_S,BC_N)

[T,Tres] = FVM_GS_ext_mesh(Nx+1,Ny+1,1.0,1e4,1e-8,T);
                   

figure('Name','Convergence Plot for Scaled Residuals',...
       'Position',[100 100 500 300])

nlist = 10:n;
semilogy(nlist,ures(nlist),'-b',nlist,vres(nlist),'-r',...
         nlist,pres(nlist),'-g')
legend('u residual','v residual','p residual')
xlabel('Iteration')
ylabel('Scaled Residual')

FVM_Vplot(Nx,Ny,x,xu,y,H,u(iu,Ju),v(Iv,jv),p(Ip,Jp),U)

FVM_Tplot(Nx,Ny,x,y,L,H,rho,cp,u(iu,Ju),T(Ip,Jp),U,Ti,Tw,qw,BC_N)