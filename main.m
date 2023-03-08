clear;
close;


%% Definições

% Número de bits
num_b = 10000000;
% Faixa de Eb/N0
Eb_N0_dB = 0:1:20;
% Faixa de Eb/N0 linearizada
Eb_N0_lin = 10 .^ (Eb_N0_dB/10);
% Viterbi profundidade
tb = 16;
% Configuração do QAM 
M = 16;
k = log2(M);


%% Fonte

source = randi([0 1], 1, num_b);
    

%% Codificador

%%% Código convolucional 1 com razão 1/2
convolutional_one_size = num_b * 2;
trellis_convolutional_one = poly2trellis(3, [7 5]);
convolutional_one_source = convenc(source, trellis_convolutional_one);

%%% Código convolucional 2 com razão 1/3
convolutional_two_size = num_b * 3;
trellis_convolutional_two = poly2trellis(3, [4 6 7]);
convolutional_two_source = convenc(source, trellis_convolutional_two);

%%% Código convolucional GSM Full Rate
convolutional_gsm_size = num_b * 2;
trellis_convolutional_gsm = poly2trellis(5, [23 27]);
convolutional_gsm_source = convenc(source, trellis_convolutional_gsm);


%% Modulação

%%% BPSK
% Sem codificação
bpsk_mod = complex(2*source-1, 0);
% Convolucional 1
bpsk_mod_convolutional_one = complex(2*convolutional_one_source-1, 0);
% Convolucional 2
bpsk_mod_convolutional_two = complex(2*convolutional_two_source-1, 0);
% Convolucional GSM
bpsk_mod_convolutional_gsm = complex(2*convolutional_gsm_source-1, 0);

%%% 16QAM
% Sem codificação
grouped_source = reshape(source, k, []).';
decimal_source = bi2de(grouped_source)';
qam_mod = qammod(decimal_source, M);
% Convolucional 1
grouped_source = reshape(convolutional_one_source, k, []).';
decimal_source = bi2de(grouped_source)';
qam_mod_convolutional_one = qammod(decimal_source, M);
% Convolucional 2
grouped_source = reshape(convolutional_two_source, k, []).';
decimal_source = bi2de(grouped_source)';
qam_mod_convolutional_two = qammod(decimal_source, M);
% Convolucional GSM
grouped_source = reshape(convolutional_gsm_source, k, []).';
decimal_source = bi2de(grouped_source)';
qam_mod_convolutional_gsm = qammod(decimal_source, M);


%% Receptor com BSPK

%%% Sem codificação
Eb = 1; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_without_code = zeros(size(Eb_N0_lin)); 
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i) * complex(randn(1, num_b), randn(1, num_b)) * sqrt(0.5); 
   
    % Parte real do vetor recebido com ruído
    bpsk_mod_with_noise = real(bpsk_mod + noise); 
    
    % Demodulação    
    demod = zeros(1, length(bpsk_mod_with_noise));
    for x = 1:length(bpsk_mod_with_noise)
        if bpsk_mod_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    disp('BPSK sem codificação');
    disp(sum(source ~= demod));
                    
    % Contagem de erros e cálculo do BER
    ber_bpsk_without_code(i) = sum(source ~= demod) / num_b; 
end

%%% Convolucional 1
Eb = 2; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_convolutional_one = zeros(size(Eb_N0_lin)); 
mod_size = length(bpsk_mod_convolutional_one);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i) * complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
   
    % Parte real do vetor recebido com ruído
    bpsk_mod_with_noise = real(bpsk_mod_convolutional_one + noise);
    
    % Demodulação    
    demod = zeros(1, length(bpsk_mod_with_noise));
    for x = 1:length(bpsk_mod_with_noise)
        if bpsk_mod_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod,trellis_convolutional_one,tb,'trunc','hard');
                
    disp('BPSK com o código convolucional 1');
    disp(sum(source ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_one(i) = sum(source ~= decoded_sequence) /  num_b; 
end

%%% Convolucional 2
Eb = 3; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP); 
ber_bpsk_convolutional_two = zeros(size(Eb_N0_lin)); 
mod_size = length(bpsk_mod_convolutional_two);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i) * complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
 
    % Parte real do vetor recebido com ruído
    bpsk_mod_with_noise = real(bpsk_mod_convolutional_two + noise);
    
    % Demodulação    
    demod = zeros(1, length(bpsk_mod_with_noise));
    for x = 1:length(bpsk_mod_with_noise)
        if bpsk_mod_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod,trellis_convolutional_two,tb,'trunc','hard');
            
    disp('BPSK com o código convolucional 2');
    disp(sum(source ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_two(i) = sum(source ~= decoded_sequence) / num_b; 
end

%%% Convolucional GSM
Eb = 2; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_bpsk_convolutional_gsm = zeros(size(Eb_N0_lin)); 
mod_size = length(bpsk_mod_convolutional_gsm);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    n = NA(i) * complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
 
    % Parte real do vetor recebido com ruído
    bpsk_mod_with_noise = real(bpsk_mod_convolutional_gsm + n); 
    
    % Demodulação    
    demod = zeros(1, length(bpsk_mod_with_noise));
    for x = 1:length(bpsk_mod_with_noise)
        if bpsk_mod_with_noise(x) > 0
            demod(x) = 1;
        else
            demod(x) = 0;
        end
    end
    
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod,trellis_convolutional_gsm,tb,'trunc','hard');
    
    disp('BPSK com o código convolucional GSM');
    disp(sum(source ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_bpsk_convolutional_gsm(i) = sum(source ~= decoded_sequence) / num_b; 
end


%% Receptor com 4QAM

%%% Sem codificação
Eb = 2.5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_without_code = zeros(size(Eb_N0_lin)); 
mod_size = length(qam_mod);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i)*complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
   
    % Vetor recebido com ruído
    qam_mod_with_noise = qam_mod + noise;
     
    % Demodulação            
    qam_demod_dec = qamdemod(qam_mod_with_noise, M);
    qam_demod_bin = de2bi(qam_demod_dec, k)';
    demod = reshape(qam_demod_bin, size(qam_demod_bin, 1) * size(qam_demod_bin, 2), 1)';
    
    disp('QAM sem codificação');
    disp(sum(source ~= demod));
    
    % Contagem de erros e cálculo do BER
    ber_qam_without_code(i) = sum(source ~= demod) / num_b; 
end

%%% Convolucional 1
Eb = 5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_one = zeros(size(Eb_N0_lin)); 
mod_size = length(qam_mod_convolutional_one);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i)*complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
   
    % Vetores recebido
    qam_mod_with_noise = qam_mod_convolutional_one + noise;
     
    % Demodulação            
    qam_demod_dec = qamdemod(qam_mod_with_noise, M);
    qam_demod_bin = de2bi(qam_demod_dec, k)';
    demod = reshape(qam_demod_bin, size(qam_demod_bin, 1) * size(qam_demod_bin, 2), 1)';
    
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod, trellis_convolutional_one, tb, 'trunc', 'hard');
       
    disp('QAM com o código convolucional 1');
    disp(sum(source ~= decoded_sequence));   
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_one(i) = sum(source ~= decoded_sequence) / num_b; 
end

%%% Convolucional 2
Eb = 7.5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_two = zeros(size(Eb_N0_lin));
mod_size = length(qam_mod_convolutional_two);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i) * complex(randn(1, mod_size), randn(1, mod_size)) * sqrt(0.5); 
   
    % Vetores recebido
    qam_mod_with_noise = qam_mod_convolutional_two + noise;
     
    % Demodulação            
    qam_demod_dec = qamdemod(qam_mod_with_noise, M);
    qam_demod_bin = de2bi(qam_demod_dec, k)';
    demod = reshape(qam_demod_bin, size(qam_demod_bin, 1) * size(qam_demod_bin, 2), 1)';
        
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod, trellis_convolutional_two, tb, 'trunc', 'hard');
       
    disp('QAM com o código convolucional 2');
    disp(sum(source ~= decoded_sequence));  
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_two(i) = sum(source ~= decoded_sequence) / num_b; 
end

%%% Convolucional GSM
Eb = 5; 
NP = Eb ./ (Eb_N0_lin); 
NA = sqrt(NP);
ber_qam_convolutional_gsm = zeros(size(Eb_N0_lin)); 
mod_size = length(qam_mod_convolutional_gsm);
for i = 1:length(Eb_N0_lin)
    % Vetor de ruído complexo com desvio padrão igual a uma posição do vetor NA
    noise = NA(i)*complex(randn(1, mod_size), randn(1, mod_size))*sqrt(0.5); 
   
    % Vetores recebido
    qam_mod_with_noise = qam_mod_convolutional_gsm + noise;

    % Demodulação            
    qam_demod_dec = qamdemod(qam_mod_with_noise, M);
    qam_demod_bin = de2bi(qam_demod_dec, k)';
    demod = reshape(qam_demod_bin, size(qam_demod_bin, 1) * size(qam_demod_bin, 2), 1)';
        
    % Decodificação usando Viterbi
    decoded_sequence = vitdec(demod, trellis_convolutional_gsm, tb, 'trunc', 'hard');
        
    disp('QAM com o código convolucional GSM');
    disp(sum(source ~= decoded_sequence));
    
    % Contagem de erros e cálculo do BER
    ber_qam_convolutional_gsm(i) = sum(source ~= decoded_sequence) / num_b; 
end


%% Avaliação resultado

% BER teórico para comparação
qam_name = strcat(int2str(M), 'QAM');
semilogy(Eb_N0_dB, ber_qam_convolutional_two, 'b', ...
    Eb_N0_dB, ber_qam_convolutional_one, 'g', ...
    Eb_N0_dB, ber_qam_convolutional_gsm, 'c', ...
    Eb_N0_dB, ber_qam_without_code, 'r', ...,
    Eb_N0_dB, ber_bpsk_convolutional_two, 'b', ...,
    Eb_N0_dB, ber_bpsk_convolutional_one, 'g', ...
    Eb_N0_dB, ber_bpsk_convolutional_gsm, 'c', ...
    Eb_N0_dB, ber_bpsk_without_code, 'r', ...
    'LineWidth', 2, 'MarkerSize', 10);
xlabel('Eb/N0 (dB)');
ylabel('BER');
legend(strcat(qam_name, ' com código convolucional de razão 1/3'), ...
    strcat(qam_name, ' com código convolucional de razão 1/2'), ...
    strcat(qam_name, ' GSM'), ...
    strcat(qam_name, ' sem codificação'), ...
    'BPSK com código convolucional de razão 1/3', ...
    'BPSK com código convolucional de razão 1/2', ...
    'BPSK com código convolucional do GSM', ...
    'BPSK sem codificação');
