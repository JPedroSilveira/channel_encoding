clear;
close;

%% Definições
% Número de bits
num_b = 100000;
% Faixa de Eb/N0
Eb_N0_dB = 0:1:9;
% Faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);

%% Fonte
info = randi(2, 1, num_b) - 1;
% Adiciona um sufixo de zeros para finalizar as máquinas convolucionais no
% estado inicial
info(end-9:end) = zeros(1, 10);
    
%% Codificador

%%% Código Convolucional 1 com razão 1/2 não recursivo e não sistemático
convolutional_machine_one = containers.Map();

% node = { edge_0 (next node, first output, second output), edge_1 (...) }
convolutional_machine_one('00') = {{'00', 0, 0}, {'10', 1, 1}};
convolutional_machine_one('10') = {{'01', 1, 0}, {'11', 0, 1}};
convolutional_machine_one('11') = {{'01', 0, 1}, {'11', 1, 0}};
convolutional_machine_one('01') = {{'00', 1, 1}, {'10', 0, 0}};

current_state = '00';
convolutional_one_size = num_b * 2;
convolutional_one_info = zeros(1, convolutional_one_size);
count = 1;

for i = 1:num_b
    state_node = convolutional_machine_one(current_state);
    state_edge = state_node{info(i) + 1};
    current_state = state_edge{1};
    convolutional_one_info(count) = state_edge{2};
    count = count + 1;
    convolutional_one_info(count) = state_edge{3};
    count = count + 1;
end

%%% Código Convolucional 2 com razão 1/3 não recursivo e não sistemático
convolutional_machine_two = containers.Map();

% node = { edge_0 (next node, first output, second output), edge_1 (...) }
convolutional_machine_two('00') = {{'00', 0, 0, 0}, {'10', 1, 1, 1}};
convolutional_machine_two('10') = {{'01', 1, 1, 1}, {'11', 0, 0, 0}};
convolutional_machine_two('11') = {{'01', 0, 0, 1}, {'11', 1, 1, 0}};
convolutional_machine_two('01') = {{'00', 1, 1, 0}, {'10', 0, 0, 1}};

current_state = '00';
convolutional_two_size = num_b * 3;
convolutional_two_info = zeros(1, convolutional_two_size);
count = 1;

for i = 1:num_b
    state_node = convolutional_machine_two(current_state);
    state_edge = state_node{info(i) + 1};
    current_state = state_edge{1};
    convolutional_two_info(count) = state_edge{2};
    count = count + 1;
    convolutional_two_info(count) = state_edge{3};
    count = count + 1;
    convolutional_two_info(count) = state_edge{4};
    count = count + 1;
end

%%% Código Convolucional GSM Full Rate
convolutional_machine_gsm = containers.Map();

% node = { edge_0 (next node, first output, second output), edge_1 (...) }
convolutional_machine_gsm('a') = {{'a', 0, 0}, {'i', 1, 1}}; 
convolutional_machine_gsm('b') = {{'a', 1, 1}, {'i', 0, 0}};
convolutional_machine_gsm('c') = {{'b', 1, 1}, {'j', 0, 0}};
convolutional_machine_gsm('d') = {{'b', 0, 0}, {'j', 1, 1}};
convolutional_machine_gsm('e') = {{'c', 0, 0}, {'k', 1, 1}};
convolutional_machine_gsm('f') = {{'c', 1, 1}, {'k', 0, 0}};
convolutional_machine_gsm('g') = {{'d', 1, 1}, {'l', 0, 0}};
convolutional_machine_gsm('h') = {{'d', 0, 0}, {'l', 1, 1}};
convolutional_machine_gsm('i') = {{'e', 0, 1}, {'m', 1, 0}};
convolutional_machine_gsm('j') = {{'e', 1, 0}, {'m', 0, 1}};
convolutional_machine_gsm('k') = {{'f', 1, 0}, {'n', 0, 1}};
convolutional_machine_gsm('l') = {{'f', 0, 1}, {'n', 1, 0}};
convolutional_machine_gsm('m') = {{'g', 0, 1}, {'o', 1, 0}};
convolutional_machine_gsm('n') = {{'g', 1, 0}, {'o', 0, 1}};
convolutional_machine_gsm('o') = {{'h', 1, 0}, {'p', 0, 1}};
convolutional_machine_gsm('p') = {{'h', 0, 1}, {'p', 1, 0}};

current_state = 'a';
convolutional_gsm_size = num_b * 2;
convolutional_gsm_info = zeros(1, convolutional_two_size);
count = 1;

for i = 1:num_b
    state_node = convolutional_machine_gsm(current_state);
    state_edge = state_node{info(i) + 1};
    current_state = state_edge{1};
    convolutional_gsm_info(count) = state_edge{2};
    count = count + 1;
    convolutional_gsm_info(count) = state_edge{3};
    count = count + 1;
end

%% Modulação

%%% BPSK
% Sem codificação
info_bpsk = complex(2*info-1, 0);
% Convolucional 1
convolutional_one_info_bpsk = complex(2*convolutional_one_info-1, 0);
% Convolucional 2
convolutional_two_info_bpsk = complex(2*convolutional_two_info-1, 0);
% Convolucional GSM
convolutional_gsm_info_bpsk = complex(2*convolutional_gsm_info-1, 0);

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

% Convolucional 1
qam_size = convolutional_one_size / 2;
convolutional_one_info_4qam_I = zeros(1, qam_size);
convolutional_one_info_4qam_Q = zeros(1, qam_size);
count = 1;
for i = 1:2:(convolutional_one_size-1)
    if convolutional_one_info(i) == 0
        if convolutional_one_info(i+1) == 0 %% 00
            convolutional_one_info_4qam_I(count) = sqrt(2)/2;
            convolutional_one_info_4qam_Q(count) = sqrt(2)/2;
        else %% 01
            convolutional_one_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_one_info_4qam_Q(count) = sqrt(2)/2;
        end
    else 
        if convolutional_one_info(i+1) == 0 %% 10
            convolutional_one_info_4qam_I(count) = sqrt(2)/2;
            convolutional_one_info_4qam_Q(count) = -1 * sqrt(2)/2;
        else %% 11
            convolutional_one_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_one_info_4qam_Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end

% Convolucional 2
qam_size = convolutional_two_size / 2;
convolutional_two_info_4qam_I = zeros(1, qam_size);
convolutional_two_info_4qam_Q = zeros(1, qam_size);
count = 1;
for i = 1:2:(convolutional_two_size-1)
    if convolutional_two_info(i) == 0
        if convolutional_two_info(i+1) == 0 %% 00
            convolutional_two_info_4qam_I(count) = sqrt(2)/2;
            convolutional_two_info_4qam_Q(count) = sqrt(2)/2;
        else %% 01
            convolutional_two_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_two_info_4qam_Q(count) = sqrt(2)/2;
        end
    else 
        if convolutional_two_info(i+1) == 0 %% 10
            convolutional_two_info_4qam_I(count) = sqrt(2)/2;
            convolutional_two_info_4qam_Q(count) = -1 * sqrt(2)/2;
        else %% 11
            convolutional_two_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_two_info_4qam_Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end

% Convolucional GSM
qam_size = convolutional_gsm_size / 2;
convolutional_gsm_info_4qam_I = zeros(1, qam_size);
convolutional_gsm_info_4qam_Q = zeros(1, qam_size);
count = 1;
for i = 1:2:(convolutional_gsm_size-1)
    if convolutional_gsm_info(i) == 0
        if convolutional_gsm_info(i+1) == 0 %% 00
            convolutional_gsm_info_4qam_I(count) = sqrt(2)/2;
            convolutional_gsm_info_4qam_Q(count) = sqrt(2)/2;
        else %% 01
            convolutional_gsm_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_gsm_info_4qam_Q(count) = sqrt(2)/2;
        end
    else 
        if convolutional_gsm_info(i+1) == 0 %% 10
            convolutional_gsm_info_4qam_I(count) = sqrt(2)/2;
            convolutional_gsm_info_4qam_Q(count) = -1 * sqrt(2)/2;
        else %% 11
            convolutional_gsm_info_4qam_I(count) = -1 * sqrt(2)/2;
            convolutional_gsm_info_4qam_Q(count) = -1 * sqrt(2)/2;
        end
    end
    count = count + 1;
end

%% Receptor com BSPK

% Energia por bit para a modulação BPSK utilizada
%%% Es = ((-1^2 + 0^2) + (1^2 0 ^2)) / 2
%%% Es = 1
%%% r = k / n
%%% r = 1
%%% Eb = Es / (M * r)
%%% Eb = 1 / (1 * 1)
%%% Eb = 1;
Eb = 1; 
% Vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 
% Vetor de amplitudes do ruído
NA = sqrt(NP); 

%%% Sem codificação
% Pré-alocação do vetor BER
ber_bpsk_without_code = zeros(size(Eb_N0_lin)); 
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
                    
    % Contagem de erros e cálculo do BER
    ber_bpsk_without_code(i) = sum(info ~= demod) / num_b; 
end

%%% Convolucional 1
% Pré-alocação do vetor BER
ber_bpsk_convolutional_one = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_one_size), randn(1, convolutional_one_size))*sqrt(0.5); 
   
    % Vetores recebido com ruído
    real_info_with_noise = real(convolutional_one_info_bpsk + n);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
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
    zero_one_state = viterbi_machine('01');
    one_zero_state = viterbi_machine('10');
    one_one_state = viterbi_machine('11');

    decoded_sequence = zero_zero_state{1};
    
    % Remove o sufixo do código convolucional
    decoded_sequence = decoded_sequence(1:end-3);
                
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_one(i) = sum(info ~= decoded_sequence) /  num_b; 
end

%%% Convolucional 2
% Pré-alocação do vetor BER
ber_bpsk_convolutional_two = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_two_size), randn(1, convolutional_two_size))*sqrt(0.5); 
 
    % Recupera a informação com ruído (sinal da parte real)
    real_info_with_noise = real(convolutional_two_info_bpsk + n);
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:3:length(demod)-2
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);
        third_bit = demod(x + 2);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = 4; % Valor máximo possível é 3
        zero_one_difference = 4;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 0
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
            if third_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 4 || zero_one_difference ~= 4
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
        one_zero_difference = 4;
        one_one_difference = 4;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 4 || one_one_difference ~= 4
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
        zero_zero_difference = 4;
        zero_one_difference = 4;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 1
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
            if third_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 4 || zero_one_difference ~= 4
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
        one_zero_difference = 4;
        one_one_difference = 4;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 4 || one_one_difference ~= 4
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
    decoded_sequence = decoded_sequence(1:end-4);
                
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_two(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% TODO: Receptor BPSK convolucional GSM
%%% Convolucional 2
% Pré-alocação do vetor BER
ber_bpsk_convolutional_gsm = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_gsm_size), randn(1, convolutional_gsm_size))*sqrt(0.5); 
 
    % Recupera a informação com ruído (sinal da parte real)
    real_info_with_noise = real(convolutional_gsm_info_bpsk ); % + n
    
    %%% Demodulação    
    demod = zeros(1, length(real_info_with_noise));
    for x = 1:length(real_info_with_noise)
        if real_info_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('a') = {[], 0};
    viterbi_machine('b') = {};
    viterbi_machine('c') = {};
    viterbi_machine('d') = {};
    viterbi_machine('e') = {};
    viterbi_machine('f') = {};
    viterbi_machine('g') = {};
    viterbi_machine('h') = {};
    viterbi_machine('i') = {};
    viterbi_machine('j') = {};
    viterbi_machine('k') = {};
    viterbi_machine('l') = {};
    viterbi_machine('m') = {};
    viterbi_machine('n') = {};
    viterbi_machine('o') = {};
    viterbi_machine('p') = {};
    
    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);

        % Armazena estados atuais
        previous_viterbi_machine = viterbi_machine;
        
        % Análise de deslocamento para a
        for a = 1
            a_difference = Inf;
            b_difference = Inf;
            a_state = previous_viterbi_machine('a');
            b_state = previous_viterbi_machine('b');
            if ~isempty(a_state)
                a_difference = a_state{2};
                if first_bit ~= 0
                    a_difference = a_difference + 1;
                end
                if second_bit ~= 0
                    a_difference = a_difference + 1;
                end
            end
            if ~isempty(b_state)
                b_difference = b_state{2};
                if first_bit ~= 1
                    b_difference = b_difference + 1;
                end
                if second_bit ~= 1
                    b_difference = b_difference + 1;
                end
            end
            if a_difference ~= Inf || b_difference ~= Inf
                if a_difference <= b_difference
                    new_bits = [a_state{1}, 0];
                    new_difference = a_difference;
                else
                    new_bits = [b_state{1}, 0];
                    new_difference = b_difference;
                end
                viterbi_machine('a') = { new_bits, new_difference };
            else
                viterbi_machine('a') = {};
            end
        end
        
        % Análise de deslocamento para b
        for b = 1
            c_difference = Inf;
            d_difference = Inf;
            c_state = previous_viterbi_machine('c');
            d_state = previous_viterbi_machine('d');
            if ~isempty(c_state)
                c_difference = c_state{2};
                if first_bit ~= 1
                    c_difference = c_difference + 1;
                end
                if second_bit ~= 1
                    c_difference = c_difference + 1;
                end
            end
            if ~isempty(d_state)
                d_difference = d_state{2};
                if first_bit ~= 0
                    d_difference = d_difference + 1;
                end
                if second_bit ~= 0
                    d_difference = d_difference + 1;
                end
            end
            if c_difference ~= Inf || d_difference ~= Inf
                if c_difference <= d_difference
                    new_bits = [c_state{1}, 0];
                    new_difference = c_difference;
                else
                    new_bits = [d_state{1}, 0];
                    new_difference = d_difference;
                end
                viterbi_machine('b') = { new_bits, new_difference };
            else
                viterbi_machine('b') = {};
            end
        end
        
        % Análise de deslocamento para c
        for c = 1
            e_difference = Inf;
            f_difference = Inf;
            e_state = previous_viterbi_machine('e');
            f_state = previous_viterbi_machine('f');
            if ~isempty(e_state)
                e_difference = e_state{2};
                if first_bit ~= 0
                    e_difference = e_difference + 1;
                end
                if second_bit ~= 0
                    e_difference = e_difference + 1;
                end
            end
            if ~isempty(f_state)
                f_difference = f_state{2};
                if first_bit ~= 1
                    f_difference = f_difference + 1;
                end
                if second_bit ~= 1
                    f_difference = f_difference + 1;
                end
            end
            if e_difference ~= Inf || f_difference ~= Inf
                if e_difference <= f_difference
                    new_bits = [e_state{1}, 0];
                    new_difference = e_difference;
                else
                    new_bits = [f_state{1}, 0];
                    new_difference = f_difference;
                end
                viterbi_machine('c') = { new_bits, new_difference };
            else
                viterbi_machine('c') = {};
            end
        end
        
    end
    
    % Seleciona o valor atual do estado 00
    zero_zero_state = viterbi_machine('00');
    decoded_sequence = zero_zero_state{1};
    
    % Remove o sufixo do código convolucional
    decoded_sequence = decoded_sequence(1:end-4);
                
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_two(i) = sum(info ~= decoded_sequence) / num_b; 
end

%% Receptor com 4-QAM
% Energia por bit para a modulação 4-QAM utilizada
%%% (2^(1/2)/)^2 = 0.5
%%% (-(2^(1/2))/2)^2 = 0.5
%%% Es = (0.5 + 0.5) + (0.5 + 0.5) + (0.5 + 0.5) + (0.5 + 0.5)
%%% Es = 4
%%% r = k / n
%%% r = 2
%%% Eb = Es / (M * r)
%%% Eb = 4 / (1 * 2)
%%% Eb = 2;
Eb = 2; 
% Vetor de potências do ruído
NP = Eb ./ (Eb_N0_lin); 
% Vetor de amplitudes do ruído
NA = sqrt(NP); 

%%% Sem codificação
% Pré-alocação do vetor BER
ber_4qam_without_code = zeros(size(Eb_N0_lin)); 
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
                
    % Contagem de erros e cálculo do BER
    ber_4qam_without_code(i) = sum(info ~= demod) / num_b; 
end

%%% Convolucional 1
% Pré-alocação do vetor BER
ber_4qam_convolutional_one = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_one_size / 2), randn(1, convolutional_one_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(convolutional_one_info_4qam_I + n); 
    Q_with_noise = real(convolutional_one_info_4qam_Q + n); 
     
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
    
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:2:length(demod)-1
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
                
    % Contagem de erros e cálculo do BER
    ber_4qam_convolutional_one(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% Convolucional 2
% Pré-alocação do vetor BER
ber_4qam_convolutional_two = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i)*complex(randn(1, convolutional_two_size / 2), randn(1, convolutional_two_size / 2))*sqrt(0.5); 
   
    % Vetores recebido
    I_with_noise = real(convolutional_two_info_4qam_I + n); 
    Q_with_noise = real(convolutional_two_info_4qam_Q + n); 
     
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
        
    %%% Decodificação usando Viterbi
    % Sequências atuais e distância total para cada estado
    viterbi_machine = containers.Map();
    viterbi_machine('00') = {[], 0};
    viterbi_machine('10') = {};
    viterbi_machine('11') = {};
    viterbi_machine('01') = {};

    % Passa pelo input e determina as sequências possíveis para a saída
    for x = 1:3:length(demod)-2
        % Bits sendo lidos
        first_bit = demod(x);
        second_bit = demod(x + 1);
        third_bit = demod(x + 2);

        % Armazena estados atuais
        viterbi_machine_00 = viterbi_machine('00');
        viterbi_machine_01 = viterbi_machine('01');
        viterbi_machine_10 = viterbi_machine('10');
        viterbi_machine_11 = viterbi_machine('11');
                
        % Análise de deslocamento para 00
        zero_zero_difference = 4; % Valor máximo possível é 3
        zero_one_difference = 4;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 0
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 0
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
            if third_bit ~= 0
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 4 || zero_one_difference ~= 4
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
        one_zero_difference = 4;
        one_one_difference = 4;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 1
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 4 || one_one_difference ~= 4
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
        zero_zero_difference = 4;
        zero_one_difference = 4;
        if ~isempty(viterbi_machine_00)
            zero_zero_difference = 0;
            if first_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if second_bit ~= 1
                zero_zero_difference = zero_zero_difference + 1;
            end
            if third_bit ~= 1
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
            if third_bit ~= 1
                zero_one_difference = zero_one_difference + 1;
            end
        end
        if zero_zero_difference ~= 4 || zero_one_difference ~= 4
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
        one_zero_difference = 4;
        one_one_difference = 4;
        if ~isempty(viterbi_machine_10)
            one_zero_difference = 0;
            if first_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if second_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
            if third_bit ~= 0
                one_zero_difference = one_zero_difference + 1;
            end
        end
        if ~isempty(viterbi_machine_11)
            one_one_difference = 0;
            if first_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if second_bit ~= 1
                one_one_difference = one_one_difference + 1;
            end
            if third_bit ~= 0
                one_one_difference = one_one_difference + 1;
            end
        end
        if one_zero_difference ~= 4 || one_one_difference ~= 4
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
    decoded_sequence = decoded_sequence(1:end-4);
                
    % Contagem de erros e cálculo do BER
    ber_4qam_convolutional_two(i) = sum(info ~= decoded_sequence) / num_b; 
end

%%% TODO: Receptor 4QAM convolucional GSM

%% Avaliação resultado

% BER teórico para comparação
ber_theoretical = 0.5*erfc(sqrt(2*Eb_N0_lin)/sqrt(2)); 

semilogy(Eb_N0_dB, ber_4qam_without_code, 'x', Eb_N0_dB, ber_4qam_convolutional_one, 'x', Eb_N0_dB, ber_4qam_convolutional_two, 'x', Eb_N0_dB, ber_bpsk_without_code, 'x', Eb_N0_dB, ber_bpsk_convolutional_one, 'x', Eb_N0_dB, ber_bpsk_convolutional_two, 'x', Eb_N0_dB, ber_theoretical, 'r', 'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend('4QAM sem codificação', '4QAM com código convolucional de razão 1/2', '4QAM com código convolucional de razão 1/3', 'BPSK sem codificação', 'BPSK com código convolucional de razão 1/2', 'BPSK com código convolucional de razão 1/3', 'Teórico');
