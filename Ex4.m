%% Lista B - Ex 4
clc,clear all,close all

%% Dados 

A = [ 0  1  0  0;
     -1  1  0  1;
     -1 -1 -1  1;
      1  0 -1  0];
B1 = [0  1  0;
      1 -1  0;
      0 -1 -1;
      1  1  0];
B2 = [ 2 -2;
      -2  0;
      -2  0;
       0  2];
C1 = [0  -1  2  1];
C2 = [-1  0 -1  0;
       1  1  2 -2];
D11 = [0 0 0];
D12 = [0 1];
D21 = [0  1  1;
       1  1 -1];

% Configurações
n = size(A,1);
r = size(B1',1);
m = size(B2',1);
p = size(C1,1);
q = size(C2,1);
output.cpusec_m = clock;
LMI = [];
LMI2 = [];

%% H2
display('a)')
% Criação das Variaveis
P = sdpvar(n,n,'symmetric');
H = sdpvar(n,n,'symmetric');
Z = sdpvar(p,p,'symmetric');
X = sdpvar(n,n,'full');
Y = sdpvar(n,n,'full');
L = sdpvar(m,n,'full');
F = sdpvar(n,q,'full');
Q = sdpvar(n,n,'full');
R = sdpvar(m,q,'full'); 
S = sdpvar(n,n,'full');
J = sdpvar(n,n,'full');
mu = sdpvar(1,1);

% Criação das LMIs
L11 = P;
L12 = J;
L13 = A*X + B2*L;
L14 = A + B2*R*C2;
L15 = B1 + B2*R*D21;
L22 = H;
L23 = Q;
L24 = Y*A + F*C2;
L25 = Y*B1 + F*D21;
L33 = X + X' - P;
L34 = eye(n) + S' - J;
L35 = zeros(n,r);
L44 = Y + Y' - H;
L45 = zeros(n,r);
L55 = eye(r);
L1 = [L11  L12  L13  L14  L15;
      L12' L22  L23  L24  L25;
      L13' L23' L33  L34  L35;
      L14' L24' L34' L44  L45;
      L15' L25' L35' L45' L55];
LMI = [LMI, L1 >= 0];

V11 = Z;
V12 = C1*X + D12*L;
V13 = C1 + D12*R*C2;
V22 = X + X' - P;
V23 = eye(n) + S' - J;
V33 = Y + Y' - H;
V1 = [V11  V12  V13;
      V12' V22  V23;
      V13' V23' V33];
LMI = [LMI, V1 >= 0];
LMI = [LMI, mu > trace(Z)];
LMI = [LMI, D11 + D12*R*D21 == 0];

% Solução a LMIs
sol_h2 = solvesdp(LMI,mu,sdpsettings('verbose',0,'solver','sedumi'));

% Resíduo
res_h2 = min(checkset(LMI));

% V e U, tomando U = eye(4);
U = eye(4);
V = value(S) - value(Y)*value(X);

% Controladores
X1 = [inv(value(V)) (-inv(value(V))*value(Y)*B2);
      zeros(m,n) eye(m)];
X2 = [(value(Q)-value(Y)*A*value(X)) value(F);
      value(L) value(R)];
X3 = [    inv(value(U))         zeros(n,q);
      -C2*value(X)*inv(value(U)) eye(q)];
cont_h2 = X1*X2*X3;

% Ganho estabilizante de realimentação de saída
display('A norma H2 é:')
sqrt(double(mu))

display('O controlador é:')

Ac = [cont_h2(1,1) cont_h2(1,2) cont_h2(1,3) cont_h2(1,4);
      cont_h2(2,1) cont_h2(2,2) cont_h2(2,3) cont_h2(2,4);
      cont_h2(3,1) cont_h2(3,2) cont_h2(3,3) cont_h2(3,4);
      cont_h2(4,1) cont_h2(4,2) cont_h2(4,3) cont_h2(4,4)]
Bc = [cont_h2(1,5) cont_h2(1,6);
      cont_h2(2,5) cont_h2(2,6);
      cont_h2(3,5) cont_h2(3,6);
      cont_h2(4,5) cont_h2(4,6)]
Cc = [cont_h2(5,1) cont_h2(5,2) cont_h2(5,3) cont_h2(5,4);
      cont_h2(6,1) cont_h2(6,2) cont_h2(6,3) cont_h2(6,4)]
Dd = [cont_h2(5,5) cont_h2(5,6);
      cont_h2(6,5) cont_h2(6,6)]

%% Hinf
display('b)')
% Criação das Variaveis
P = sdpvar(n,n,'symmetric');
H = sdpvar(n,n,'symmetric');
X = sdpvar(n,n,'full');
Y = sdpvar(n,n,'full');
L = sdpvar(m,n,'full');
F = sdpvar(n,q,'full');
Q = sdpvar(n,n,'full');
R = sdpvar(m,q,'full');
S = sdpvar(n,n,'full');
J = sdpvar(n,n,'full');
gamma = sdpvar(1,1);

% Criação das LMIs
T11 = P;
T12 = J;
T13 = A*X + B2*L;
T14 = A + B2*R*C2;
T15 = B1 + B2*R*D21;
T16 = zeros(n,p);
T22 = H;
T23 = Q;
T24 = Y*A + F*C2;
T25 = Y*B1 + F*D21;
T26 = zeros(n,p);
T33 = X + X' - P;
T34 = eye(n) + S' - J;
T35 = zeros(n,r);
T36 = X'*C1' + L'*D12';
T44 = Y + Y' - H;
T45 = zeros(n,r);
T46 = C1' + C2'*R'*D12';
T55 = eye(r);
T56 = D11' + D21'*R'*D12';
T66 = gamma*eye(p);
T1 = [T11  T12  T13  T14  T15  T16;
      T12' T22  T23  T24  T25  T26;
      T13' T23' T33  T34  T35  T36;
      T14' T24' T34' T44  T45  T46;
      T15' T25' T35' T45' T55  T56;
      T16' T26' T36' T46' T56' T66];
LMI2 = [LMI2,T1 >= 0];

% Solução a LMIs
sol_hinf = solvesdp(LMI2,gamma,sdpsettings('verbose',0,'solver','sedumi'));

% Resíduo
res_hinf = min(checkset(LMI2));

% V e U, tomando U = eye(4);
U = eye(4);
V = value(S) - value(Y)*value(X);

% Controladores
X1 = [inv(value(V)) (-inv(value(V))*value(Y)*B2);
      zeros(m,n) eye(m)];
X2 = [(value(Q)-value(Y)*A*value(X)) value(F);
      value(L) value(R)];
X3 = [    inv(value(U))         zeros(n,q);
      -C2*value(X)*inv(value(U)) eye(q)];
cont_hinf = X1*X2*X3;

% Ganho estabilizante de realimentação de saída
display('A norma Hinf é:')
sqrt(double(gamma))

display('O controlador é:')

Ac = [cont_hinf(1,1) cont_hinf(1,2) cont_hinf(1,3) cont_hinf(1,4);
      cont_hinf(2,1) cont_hinf(2,2) cont_hinf(2,3) cont_hinf(2,4);
      cont_hinf(3,1) cont_hinf(3,2) cont_hinf(3,3) cont_hinf(3,4);
      cont_hinf(4,1) cont_hinf(4,2) cont_hinf(4,3) cont_hinf(4,4)]
Bc = [cont_hinf(1,5) cont_hinf(1,6);
      cont_hinf(2,5) cont_hinf(2,6);
      cont_hinf(3,5) cont_hinf(3,6);
      cont_hinf(4,5) cont_hinf(4,6)]
Cc = [cont_hinf(5,1) cont_hinf(5,2) cont_hinf(5,3) cont_hinf(5,4);
      cont_hinf(6,1) cont_hinf(6,2) cont_hinf(6,3) cont_hinf(6,4)]
Dd = [cont_hinf(5,5) cont_hinf(5,6);
      cont_hinf(6,5) cont_hinf(6,6)]