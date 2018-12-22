%% Lista B - Ex 3
clc,clear all,close all

%% Dados 

a{1} = [1 1; 5 -2];
a{2} = [-4 -7  2; 0 -4 -2; -3  2  0];
a{3} = [-1  2  7; 3 -2  0; -3 -4 -1];

b{1} = [3; 2];
b{2} = [1; 0; 1];
b{3} = [-1 -1; -1 -2; 2 -1];

c{1} = [1 0];
c{2} = [1 0 0];
c{3} = [1 0 0];

escolha = 3; % a) = 1 b) = 2 c) = 3

A = a{escolha};
B = b{escolha};
C = c{escolha};

% Configurações
n = size(A,1);
m = size(B',1);
p = size(C,1);
output.cpusec_m = clock;
LMI = [];

% Criação das Variaveis
W11  = sdpvar(p,p,'symmetric');
W22  = sdpvar((n-p),(n-p),'symmetric');
Z11  = sdpvar(m,1,'full');

k = 1;
T = k*eye(n);
T(1,:) = C;

A_ = T*A*inv(T);
B_ = T*B;

% Criação das LMIs
W0 = [W11 zeros(size(W11,1),size(W22',1));zeros(size(W22,1),size(W11',1)) W22];
LMI = [LMI, W0 >= 0];
Z0 = [Z11 zeros(size(Z11,1),(size(B_,1)-size(Z11',1)))];
L1 = A_*W0 + W0*A_' + B_*Z0 + Z0'*B_';
LMI = [LMI, L1 <= 0];

% Número de Variáveis
output.V = size(getvariables(LMI),2);


% Solução a LMIs
sol = solvesdp(LMI,[],sdpsettings('verbose',0,'solver','sedumi'));

% Tempo para solução
output.cpusec_s = sol.solvertime;

% Resíduo
res = min(checkset(LMI));

% Ganho estabilizante de realimentação de saída
display('O Ganho estabilizante de realimentação de saída é:')
L = value(Z11)*inv(value(W11))
novo_A = A + B*L*C;
display('Autovalores antes do ganho:')
eig(A)
display('Autovalores depois do ganho:')
eig(novo_A)

