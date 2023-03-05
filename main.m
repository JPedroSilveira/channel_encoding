clear;
close;

%% Definições
% Número de bits
num_b = 6;
% Faixa de Eb/N0
Eb_N0_dB = 0:1:9;
% Faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);

%% Fonte
info = randi(2, 1, num_b) - 1;
    
%% Codificador

%%% Viterbi
convolutional_machine = containers.Map();

% node = { edge_0 (next node, first output, second output), edge_1 (...) }
convolutional_machine('00') = {{'00', 0, 0}, {'10', 1, 1}};
convolutional_machine('10') = {{'01', 1, 0}, {'11', 0, 1}};
convolutional_machine('11') = {{'01', 0, 1}, {'11', 1, 0}};
convolutional_machine('01') = {{'00', 1, 1}, {'10', 0, 0}};

% Sufixo para garantir que o estado final será 00
convolutional_end = containers.Map();
convolutional_end('00') = [0, 0, 0, 0, 0, 0];
convolutional_end('10') = [0, 1, 0, 1, 1, 1];
convolutional_end('11') = [0, 1, 1, 1, 0, 0];
convolutional_end('01') = [1, 1, 0, 0, 0, 0];

current_state = '00';
viterbi_info_size = num_b * 2 + 6;
viterbi_info = zeros(1, viterbi_info_size);
count = 1;

for i = 1:num_b
    state_node = convolutional_machine(current_state);
    state_edge = state_node{info(i) + 1};
    current_state = state_edge{1};
    viterbi_info(count) = state_edge{2};
    count = count + 1;
    viterbi_info(count) = state_edge{3};
    count = count + 1;
end

viterbi_suffix = convolutional_end(current_state);
for i = 1:length(viterbi_suffix)
    viterbi_info(count) = viterbi_suffix(i);
    count = count + 1;
end

%% Modulação

%%% BPSK
% Sem codificação
info_bpsk = complex(2*info-1, 0);
% Viterbi
viterbi_info_bpsk = complex(2*viterbi_info-1, 0);

%%% 4-QAM
% Bits | Q          | I
% 00   | sqrt(2)/2  | sqrt(2)/2
% 01   | sqrt(2)/2  | -sqrt(2)/2
% 10   | -sqrt(2)/2 | sqrt(2)/2
% 11   | -sqrt(2)/2 | -sqrt(2)/2

% Sem codificação
qam_size = num_b / 2;
info_4qam_I = zeros(1, qam_size);
info_4qam_Q = zeros(1, qam_size);
count = 1;
for i = 1:2:(length(info)-1)
    if info(i) == 0
        if info(i+1) == 0 %% 00
            info_4qam_I(count) = sqrt(2)/2;
            info_4qam_Q(count) = sqrt(2)/2;
        else %% 01
            info_4qam_I(count) = -1 * sqrt(2)/2;
            info_4qam_Q(count) = sqrt(2)/2;
        end
    else 
        if info(i+1) == 0 %% 10
            info_4qam_I(count) = sqrt(2)/2;
            info_4qam_Q(count) = -1 * sqrt(2)/2;
        else %% 11
            info_4qam_I(count) = -1 * sqrt(2)/2;
            info_4qam_Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end

% Viterbi
qam_size = viterbi_info_size / 2;
viterbi_info_4qam_I = zeros(1, qam_size);
viterbi_info_4qam_Q = zeros(1, qam_size);
count = 1;
for i = 1:2:(viterbi_info_size-1)
    if viterbi_info(i) == 0
        if viterbi_info(i+1) == 0 %% 00
            viterbi_info_4qam_I(count) = sqrt(2)/2;
            viterbi_info_4qam_Q(count) = sqrt(2)/2;
        else %% 01
            viterbi_info_4qam_I(count) = -1 * sqrt(2)/2;
            viterbi_info_4qam_Q(count) = sqrt(2)/2;
        end
    else 
        if viterbi_info(i+1) == 0 %% 10
            viterbi_info_4qam_I(count) = sqrt(2)/2;
            viterbi_info_4qam_Q(count) = -1 * sqrt(2)/2;
        else %% 11
            viterbi_info_4qam_I(count) = -1 * sqrt(2)/2;
            viterbi_info_4qam_Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end

%% Receptor com BSPK

% Energia por bit para a modulação BPSK utilizada
Eb = 1; 
% Vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 
% Vetor de amplitudes do ruído
NA = sqrt(NP); 

% Sem codificação
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b), randn(1, num_b))*sqrt(0.5); 
   
    % Vetores recebido
    info_bpsk_with_noise = info_bpsk + n; 
        
    % Recupera a informação (sinal da parte real)
    real_info_with_noise = real(info_bpsk_with_noise);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    disp('Sem codificação BPSK');
    disp(sum(info ~= demod));
                
    % Contagem de erros e cálculo do BER
    % ber(i) = sum(bits ~= demod) / viterbi_info_size; 
end

% Viterbi
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, viterbi_info_size), randn(1, viterbi_info_size))*sqrt(0.5); 
   
    % Vetores recebido
    viterbi_info_bpsk_with_noise = viterbi_info_bpsk; % + n; 
        
    % Recupera a informação (sinal da parte real)
    real_info_with_noise = real(viterbi_info_bpsk_with_noise);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação   
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = 3; % Valor máximo possível é 2
        zero_one_difference = 3;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = 0;
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 3 || zero_one_difference ~= 3
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], origin_state{2} + difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = 3;
        one_one_difference = 3;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 3 || one_one_difference ~= 3
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], origin_state{2} + difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = 3;
        zero_one_difference = 3;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = 0;
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 3 || zero_one_difference ~= 3
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], origin_state{2} + difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = 3;
        one_one_difference = 3;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 3 || one_one_difference ~= 3
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], origin_state{2} + difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
    
    % Remove o sufixo do código convolucional
    decoded_sequence = decoded_sequence(1:end-3);
    
    disp('Viterbi BPSK');
    disp(sum(info ~= decoded_sequence));
                
    % Contagem de erros e cálculo do BER
    % ber(i) = sum(bits ~= demod) / viterbi_info_size; 
end

%% Receptor com 4-QAM
% Energia por bit para a modulação 4-QAM utilizada
Eb = 1; % TODO: Calcular
% Vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 
% Vetor de amplitudes do ruído
NA = sqrt(NP); 

% Sem codificação
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b / 2), randn(1, num_b / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(info_4qam_I + n);
    Q_with_noise = real(info_4qam_Q + n);
     
    %%% Demodulação            
    demod = zeros(1, length(I_with_noise) * 2);
    count = 1;
    for x = 1:length(I_with_noise)
        if I_with_noise(x) >= 0
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
    
    disp('Sem codificação 4QAM');
    disp(sum(info ~= demod));
                
    % Contagem de erros e cálculo do BER
    % ber(i) = sum(bits ~= demod) / viterbi_info_size; 
end

% Viterbi
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, viterbi_info_size / 2), randn(1, viterbi_info_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    viterbi_I_with_noise = real(viterbi_info_4qam_I + n); 
    viterbi_Q_with_noise = real(viterbi_info_4qam_Q + n); 
     
    %%% Demodulação            
    demod = zeros(1, length(viterbi_I_with_noise) * 2);
    count = 1;
    for x = 1:length(viterbi_I_with_noise)
        if viterbi_I_with_noise(x) >= 0
            if viterbi_Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 0;
            else
                demod(count) = 1;
                demod(count + 1) = 0;
            end
        else
            if viterbi_Q_with_noise(x) >= 0
                demod(count) = 0;
                demod(count + 1) = 1;
            else
                demod(count) = 1;
                demod(count + 1) = 1;
            end
        end
        count = count + 2;
    end
    
    %%% Decodificação   
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = 3; % Valor máximo possível é 2
        zero_one_difference = 3;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = 0;
            if first_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 3 || zero_one_difference ~= 3
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('00') = { [origin_state{1}, 0], origin_state{2} + difference };
        else
            viterbi_machine('00') = {};
        end
        
        % Análise de deslocamento para 01
        one_zero_difference = 3;
        one_one_difference = 3;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 3 || one_one_difference ~= 3
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('01') = { [origin_state{1}, 0], origin_state{2} + difference };
        else
            viterbi_machine('01') = {};
        end
        
        % Análise de deslocamento para 10
        zero_zero_difference = 3;
        zero_one_difference = 3;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_01)
            zero_one_difference = 0;
            if first_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
            if second_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 3 || zero_one_difference ~= 3
            if zero_zero_difference <= zero_one_difference
                origin_state = viterbi_machine_00;
                difference = zero_zero_difference;
            else
                origin_state = viterbi_machine_01;
                difference = zero_one_difference;
            end
            viterbi_machine('10') = { [origin_state{1}, 1], origin_state{2} + difference };
        else
            viterbi_machine('10') = {};
        end
        
        % Análise de deslocamento para 11        
        one_zero_difference = 3;
        one_one_difference = 3;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 3 || one_one_difference ~= 3
            if one_zero_difference <= one_one_difference
                origin_state = viterbi_machine_10;
                difference = one_zero_difference;
            else
                origin_state = viterbi_machine_11;
                difference = one_one_difference;
            end
            viterbi_machine('11') = { [origin_state{1}, 1], origin_state{2} + difference };
        else
            viterbi_machine('11') = {};
        end
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
    
    % Remove o sufixo do código convolucional
    decoded_sequence = decoded_sequence(1:end-3);
    
    disp('Viterbi 4QAM');
    disp(sum(info ~= decoded_sequence));
                
    % Contagem de erros e cálculo do BER
    % ber(i) = sum(bits ~= demod) / viterbi_info_size; 
end

%% Avaliação resultado

%% Prof
% bits aleatórios modulados em BPSK (parte real em 1 e -1)
bits = complex(2*randi(2, 1, num_b)-3, 0); 

% faixa de Eb/N0
Eb_N0_dB = 0:1:9; 

% faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);                                

% pré-alocação do vetor de BER
ber = zeros(size(Eb_N0_lin)); 

% energia por bit para a modulação BPSK utilizada
Eb = 1; 

% vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin);

% vetor de amplitudes do ruído
NA = sqrt(NP); 


for i = 1:length(Eb_N0_lin)
    % vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, num_b), randn(1, num_b))*sqrt(0.5); 
   
    % vetor recebido
    r = bits + n; 
        
    % recupera a informação (sinal da parte real)
    demod = sign(real(r)); 
            
    % contagem de erros e cálculo do BER
    ber(i) = sum(bits ~= demod) / num_b; 
end

%BER teórico para comparação
ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); 

semilogy(Eb_N0_dB, ber, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('Simulado','Teórico');
