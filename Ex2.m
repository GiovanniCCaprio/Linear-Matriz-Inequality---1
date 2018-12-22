%% Lista B - Ex 2
clc,clear all,close all

% Dados
A = [-0.9 1; 0 0.9];
B1 = [-2; 4];
C1 = [1 0];
C2 = [-1 -3];
D11 = 0.1;
D21 = 0;
N  =  50;


%% a) Hinfinito 
display('a)')

% Configurações
n = size(A,1);
r = size(B1',1);
p = size(C1,1);
q = size(C2,1);
output.cpusec_m = clock;
LMI = [];
x = [];
xf = [];
z = [];
zf = [];
x{1} = [10;-10];
xf{1} = [0; 0];
e = 0;

% Criação das Variaveis
Z  = sdpvar(n,n,'symmetric');
X  = sdpvar(n,n,'symmetric');
F  = sdpvar(p,n,'full');
L  = sdpvar(n,q,'full');
G  = sdpvar(n,n,'full');
Df = sdpvar(p,q,'full');
gamma = sdpvar(1,1);

% Criação das LMIs
L11 = Z;
L13 = Z*A;
L15 = Z*B1;
L16 = zeros(n,p);
L26 = zeros(n,p);
L22 = X;
L23 = X*A + L*C2 + G;
L24 = X*A + L*C2;
L25 = X*B1 + L*D21;
L53 = zeros(n,p);
L54 = zeros(n,q);
L55 = eye(r); 
L63 = C1' -C2'*Df' - F'; 
L64 = C1' -C2'*Df'; 
L65 = D11' - D21'*Df';
L66 = gamma*eye(p);
L1  = [L11  L11  L13  L13  L15  L16;
       L11' L22  L23  L24  L25  L26;
       L13' L23' L11  L11  L53  L63;
       L13' L24' L11' L22  L54  L64;
       L15' L25' L54' L53' L55  L65;
       L16' L26' L63' L64' L65' L66];
LMI = [LMI, L1 >= 0];

% Número de Variáveis
output.V = size(getvariables(LMI),2);

% Solução a LMIs
sol = solvesdp(LMI,gamma,sdpsettings('verbose',0,'solver','sedumi'));

% Tempo para solução
output.cpusec_s = sol.solvertime;

% Resíduo
res_hinf = min(checkset(LMI));

% Tomando U = eye; ((XZ^(-1) + U'V = I ))
V = eye(n) - (value(X)*inv(value(Z)));

% Filtro 
Af = value(G)*inv(V*value(Z))
Bf = value(L)
Cf = value(F)*inv(V*value(Z))

% Simulação temporal(POSSÍVEL ERRO NA DINÂMICA)
for k = 1:(N)
    x{k+1}  = A*x{k}   + B1*exp(-0.25*k)*cos(10*k);
    z{k}    = C1*x{k}  + D11*exp(-0.25*k)*cos(10*k);
    xf{k+1} = Af*xf{k} + Bf*(C2*x{k} + D21*exp(-0.25*k)*cos(10*k));
    zf{k}   = Cf*xf{k}  + value(Df)*(C2*x{k} + D21*exp(-0.25*k)*cos(10*k));
    e = e + (z{k} - zf{k})^2; 
end

% Valor da Norma Hinf
display('O Valor da Norma Hinf:')
sqrt(value(gamma))

% Erro Quadrático Médio
display('O Erro Quadrático Médio é:')
e/N

%% b) H2 
display('b)')

% Configurações
n = size(A,1);
r = size(B1',1);
p = size(C1,1);
q = size(C2,1);
output.cpusec_m = clock;
LMI_h2 = [];
x = [];
xf = [];
z = [];
zf = [];
x{1} = [10;-10];
xf{1} = [0; 0];
e2 = 0;
Df = 0;

% Criação das Variaveis
Z  = sdpvar(n,n,'symmetric');
X  = sdpvar(n,n,'symmetric');
F  = sdpvar(p,n,'full');
L  = sdpvar(n,q,'full');
G  = sdpvar(n,n,'full');
M  = sdpvar(r,r,'full');
mu = sdpvar(1,1);

% Criação das LMIs
L11 = Z;
L13 = Z*B1;
L22 = X;
L23 = X*B1 + L*D21;
L33 = M;
L1  = [L11  L11  L13;
       L11  L22  L23;
       L13' L23' L33];
LMI_h2 = [LMI_h2, L1 >= 0];

L11 = Z;
L13 = A'*Z;
L14 = A'*X + C2'*L' +G'; 
L15 = C1' - F';
L22 = X;
L24 = A'*X + C2'*L';
L25 = C1';
L35 = zeros(n,r);
L45 = zeros(n,r);
L55 = eye(r);
L2  = [L11  L11  L13  L14  L15;
       L11' L22  L13  L24  L25;
       L13' L13' L11  L11  L35;
       L14' L24' L11' L22  L45;
       L15' L25' L35' L45' L55];
LMI_h2 = [LMI_h2, L2 >= 0];

LMI_h2 = [LMI_h2, mu >= trace(M)];

% Número de Variáveis
output.V = size(getvariables(LMI_h2),2);

% Solução a LMIs
sol = solvesdp(LMI_h2,mu,sdpsettings('verbose',0,'solver','sedumi'));

% Tempo para solução
output.cpusec_s = sol.solvertime;

% Resíduo
res_h2 = min(checkset(LMI_h2));

% Tomando U = eye; ((XZ^(-1) + U'V = I ))
V = eye(n) - (value(X)*inv(value(Z)));

% Filtro 
Af = value(G)*inv(V*value(Z))
Bf = value(L)
Cf = value(F)*inv(V*value(Z))

% Simulação temporal(POSSÍVEL ERRO NA DINÂMICA)
for k = 1:(N)
    x{k+1}  = A*x{k}   + B1*exp(-0.25*k)*cos(10*k);
    z{k}    = C1*x{k}  + D11*exp(-0.25*k)*cos(10*k);
    xf{k+1} = Af*xf{k} + Bf*(C2*x{k} + D21*exp(-0.25*k)*cos(10*k));
    zf{k}   = Cf*xf{k}  + value(Df)*(C2*x{k} + D21*exp(-0.25*k)*cos(10*k));
    e2 = e2 + (z{k} - zf{k})^2; 
end

% Valor da Norma Hinf
display('O Valor da Norma H2:')
sqrt(double(mu))

% Erro Quadrático Médio
display('O Erro Quadrático Médio é:')
e2/N