%% Lista B - Ex 5
clc,clear all,close all

%% Dados 

A = [ 1 -2  6;
      7  2 -5;
     -9  7  4];
B = [ 2  0;
     -3 -1;
      2  1];

% Configurações
n = size(A,1);
m = size(B',1);
output.cpusec_m = clock;
LMI_a = [];
LMI_b = [];
LMI_c = [];
LMI_d = [];
LMI_e = [];

%% a) Busca pelo ganho "cheio"
display('A)')

% Criação das Variaveis
Wa = sdpvar(n,n,'symmetric');
Za = sdpvar(m,n,'full');

% Criação das LMIs
LMI_a = [LMI_a, Wa >= 0];
La = A*Wa + Wa*A' + B*Za + Za'*B'; 
LMI_a = [LMI_a, La <= 0];

% Solução a LMIs
sol_a = solvesdp(LMI_a,[],sdpsettings('verbose',0,'solver','sedumi'));

% Resíduo
res_a = min(checkset(LMI_a))

% Ganho
display('O ganho K é:')
K_a = value(Za)*inv(value(Wa))

display('Autovalores sem o ganho K:')
eig(A)
display('Autovalores com o ganho K:')
eig(A+B*K_a)

%% b) Ganho descentralizado
display('B)')

% Criação das Variaveis
W11 = sdpvar();
W13 = sdpvar();
W22 = sdpvar();
W31 = sdpvar();
W33 = sdpvar();
Wb = [W11 0 W13; 0 W22 0; W31 0 W33];

Z11 = sdpvar();
Z22 = sdpvar();
Z13 = sdpvar();
Zb = [Z11 0 Z13; 0 Z22 0];

% Criação das LMIs
LMI_b = [LMI_b, Wb >= 0];
Lb = A*Wb + Wb*A' + B*Zb + Zb'*B'; 
LMI_b = [LMI_b, Lb <= 0];

% Solução a LMIs
sol_b = solvesdp(LMI_b,[],sdpsettings('verbose',0,'solver','sedumi'));

% Resíduo
res_b = min(checkset(LMI_b))

% Ganho
display('O ganho K é:')
K_b = value(Zb)*inv(value(Wb))

%% c) 
display('C)')

% Criação das Variaveis
Wc = sdpvar(n,n,'symmetric');
X = sdpvar(n,n,'symmetric');
Z = sdpvar(m,n,'full');

e{1} = 0.001; %inicial
e{2} = 0.01;
e{3} = 0.1;
e{4} = 1;
e{5} = 10;
e{6} = 100;

for k = 1:6
% Criação das LMIs
L11 = A*X + X'*A' + B*Z + Z'*B'; 
L12 = Wc - X' + e{k}*A*X + e{k}*B*Z;
L22 = -e{k}*(X + X');
L1 = [L11  L12;
      L12' L22];
LMI_c = [LMI_c, L1 <= 0];

% Solução a LMIs
sol_c = solvesdp(LMI_c,[],sdpsettings('verbose',0,'solver','sedumi'));

Epsilon = e{k}
% Resíduo
res_c = min(checkset(LMI_c))

% Ganho
display('O ganho Kc é:')
K_c = value(Z)*inv(value(X)) % Ganho K final

LMI_c = []; %Zera pro próximo cálculo
end

%% d) Ganho descentralizado
display('D)')


% Dados
e{1} = 0.001; %inicial
e{2} = 0.01;
e{3} = 0.1;
e{4} = 1;
e{5} = 10;
e{6} = 100;

for k = 1:6
    % Criação das Variaveis
    Wd = sdpvar(n,n,'symmetric');
    Z11 = sdpvar();
    Z22 = sdpvar();
    Z13 = sdpvar();
    Z = [Z11 0 Z13; 0 Z22 0];
    %X = sdpvar(n,n,'diagonal');
    X = [sdpvar() 0 sdpvar();0 sdpvar() 0 ;sdpvar() 0 sdpvar()]; % descentralizado.

    % Criação das LMIs
    L11 = A*X + X'*A' + B*Z + Z'*B'; 
    L12 = Wd - X' + e{k}*A*X + e{k}*B*Z;
    L22 = -e{k}*(X + X');
    L1 = [L11  L12;
          L12' L22];
    LMI_d = [LMI_d, L1 <= 0];

    % Solução a LMIs
    sol_d = solvesdp(LMI_d,[],sdpsettings('verbose',0,'solver','sedumi'));

    Epsilon = e{k}
    % Resíduo
    res_d = min(checkset(LMI_d))
    
    % Ganho
    K_d{k} = value(Z)*inv(value(X))
    
    LMI_d = []; %Zera pro próximo cálculo
end

%% e) 
display('E)')

    % 1º Estágio
    display('1º ESTÁGIO')
% Dados
e{1} = 0.001; %inicial
e{2} = 0.01;
e{3} = 0.1;
e{4} = 1;
e{5} = 10;
e{6} = 100;

for k = 1:6
    % Criação das Variaveis
    We = sdpvar(n,n,'symmetric');
    Z = sdpvar(m,n,'full');
    X = sdpvar(n,n,'symmetric');

    % Criação das LMIs
    L11 = A*X + X'*A' + B*Z + Z'*B'; 
    L12 = We - X' + e{k}*A*X + e{k}*B*Z;
    L22 = -e{k}*(X + X');
    L1 = [L11  L12;
          L12' L22];
    LMI_e = [LMI_e, L1 <= 0];

    % Solução a LMIs
    sol_e = solvesdp(LMI_e,[],sdpsettings('verbose',0,'solver','sedumi'));
    
    Epsilon = e{k}
    % Resíduo
    res_e = min(checkset(LMI_e))
    % Ganho
    K_e{k} = value(Z)*inv(value(X)); % Ganho Cheio
    
    LMI_e = []; %Zera pro próximo cálculo
end

    % 2º Estágio
    display('2º ESTÁGIO')
    
% Dados
C = eye(n);

for w = 1:6
% Criação das Variaveis
P = sdpvar(n,n,'symmetric');
F = sdpvar(n,n,'full');
G = sdpvar(n,n,'full');
H = sdpvar(m,m,'diagonal');
%J = sdpvar(m,n,'full'); 
J11 = sdpvar(1,1,'full');
J22 = sdpvar(1,1,'full');
J13 = sdpvar(1,1,'full');
J = [J11 0 J13; 0 J22 0];

% Criação das LMIs
LMI_e = [LMI_e, P >= 0];

E11 = A'*F + F'*A + K_e{w}'*B'*F' + F*B*K_e{w}; 
E12 = P - F +A'*G' + K_e{w}'*B'*G';
E13 = F*B + C'*J' - K_e{w}'*H';
E22 = -G -G';
E23 = G*B;
E33 = -H -H';
E = [E11  E12  E13;
     E12' E22  E23;
     E13' E23' E33];
 
LMI_e = [LMI_e, E <= 0];

% Solução a LMIs
sol_e = solvesdp(LMI_e,[],sdpsettings('verbose',0,'solver','sedumi'));

Epsilon = e{w}
% Resíduo
res_e = min(checkset(LMI_e))

% Ganho 
display('O ganho L é:')
L = inv(value(H))*value(J)

LMI_e = []; %Zera pro próximo cálculo
end

