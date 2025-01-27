Funcionalidades

    Condicao_de_Fase: Calcula a fase total de um sistema de controle e determina a fase necessária para correção com base em um polo desejado.
    Condicao_de_Magnitude: Calcula o valor de KK para um sistema de controle em malha fechada e a sua fase.
    Calculo_T_s: Calcula a função de transferência em malha fechada T(s)T(s) e os polos do sistema.
    Resposta_degrau: Plota a resposta ao degrau de um sistema de controle.
    Lugar_geometrico_das_raizes: Plota o lugar geométrico das raízes (Root Locus) para o sistema de controle.
    Diagrama_Bode: Plota o diagrama de Bode do sistema de controle, usando dados de um arquivo Excel ou de uma função de transferência fornecida.
    Exemplo de Uso

Para usar a classe controle_classico e suas funções, siga o exemplo abaixo:

    # Defina o polo desejado e as funções de transferência
    s = symbols('s')
    Polo_desejado = -0.25 + 1j
    num_H_s = [10]
    den_H_s = [s, s - 1, s + 1, s + 2]
    # Calcule a condição de fase
    controle_classico.Condicao_de_Fase(num_H_s, den_H_s, Polo_desejado)
