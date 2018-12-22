%% Lista B - Ex 1
clc,clear all,close all

% Dados

A{1} = [0 1 0; 0 0 1; -1 -3 -2];
A{2} = [0 1 0; 0 0 1; -1 -22 -2];
H = [1 -1;-1  1]; 

% Configurações

order = size(A{1},1);
output.cpusec_m = clock;
val = [];
gamma = 0.001;
k = 1;
p = 0;
while p >=0
    % Criação das Variaveis
    P1 = sdpvar(order,order,'symmetric'); % Lyapunov > 0
    P2 = sdpvar(order,order,'symmetric'); % Lyapunov > 0
    G1 = sdpvar(order,order,'full'); 
    G2 = sdpvar(order,order,'full');
    F = sdpvar(order,order,'full');
 
    % Variação a cada novo gamma dos h
    h1 = gamma*H(:,1);
    h2 = gamma*H(:,2);

    % Criação das LMIs
    LMIs = []; % novas LMIs a cada iteração
    LMIs = [LMIs, P1 >= 0];
    LMIs = [LMIs, P2 >= 0];
    
    % B1*alpha1^2
    L1_11 = h1(1,1)*P1 + h1(2,1)*P2 + F*A{1} + A{1}'*F'; 
    L1_12 = P1 -F + A{1}'*G1';
    L1_22 = -G1-G1';
    L1 = [L1_11 L1_12; L1_12' L1_22];
    LMIs = [LMIs, L1 <= 0];
   
    % B1*alpha1*alpha2
    L2_11 = 2*(h1(1,1)*P1 + h1(2,1)*P2) + (F*A{1} + F*A{2}) + (A{1}'*F' + A{2}'*F');
    L2_12 = (P1 + P2) -(2*F) +(A{1}'*G2' + A{2}'*G1');
    L2_22 = -G2-G1-G2'-G1';
    L2 = [L2_11 L2_12; L2_12' L2_22];
    LMIs = [LMIs, L2 <= 0];
    
    % B1*alpha2^2
    L3_11 = h1(1,1)*P1 + h1(2,1)*P2 + F*A{2} + A{2}'*F'; 
    L3_12 = P2 -F + A{2}'*G2';
    L3_22 = -G2-G2';
    L3 = [L3_11 L3_12; L3_12' L3_22];
    LMIs = [LMIs, L3 <= 0];
    
    % B2*alpha1^2
    L4_11 = h2(1,1)*P1 + h2(2,1)*P2 + F*A{1} + A{1}'*F'; 
    L4_12 = P1 -F + A{1}'*G1';
    L4_22 = -G1-G1';
    L4 = [L4_11 L4_12; L4_12' L4_22];
    LMIs = [LMIs, L4 <= 0];
    
    % B2*alpha1*alpha2
    L5_11 = 2*(h2(1,1)*P1 + h2(2,1)*P2) + (F*A{1} + F*A{2}) + (A{1}'*F' + A{2}'*F');
    L5_12 = (P1 + P2) -(2*F) +(A{1}'*G2' + A{2}'*G1');
    L5_22 = -G2-G1-G2'-G1';
    L5 = [L5_11 L5_12; L5_12' L5_22];
    LMIs = [LMIs, L5 <= 0];
    
    % B2*alpha2^2
    L6_11 = h2(1,1)*P1 + h2(2,1)*P2 + F*A{2} + A{2}'*F'; 
    L6_12 = P2 -F + A{2}'*G2';
    L6_22 = -G2-G2';
    L6 = [L6_11 L6_12; L6_12' L6_22];
    LMIs = [LMIs, L6 <= 0];

    % Número de Variáveis
    output.V = size(getvariables(LMIs),2);

    % Solução a LMIs
    sol = solvesdp(LMIs,[],sdpsettings('verbose',0,'solver','sedumi'));

    % Tempo para solução
    output.cpusec_s = sol.solvertime;

    % Resíduo
    p = min(checkset(LMIs));
    val{k} = p; % salva os residuos
    k = k + 1; % nova iteração
    gamma = 0.001 + k*0.01; % gamma novo
end

display('O gamma máximo é:')
gamma_max = gamma - 0.01