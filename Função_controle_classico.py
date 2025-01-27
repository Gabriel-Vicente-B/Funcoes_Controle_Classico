import numpy as np
from sympy import symbols, simplify, fraction, expand
import matplotlib.pyplot as plt
from scipy import signal
from control.matlab import tf, rlocus
import pandas as pd
import sympy as sp
s = symbols('s')

class controle_classico:
    def Condicao_de_Fase(num_H_s,den_H_s,polo_desejado):
        

        num_ft = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in num_H_s]
        den_ft = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in den_H_s]
        num = [complex(val) for val in num_ft]
        den= [complex(val) for val in den_ft]

        fase_total_num=0
        fase_total_den=0
        for termo_num in num:
            if termo_num == int: 
                fase_total_num=0
                pass
            fase_total_num += np.degrees(np.angle(termo_num))
        for termo_den in den:
            if termo_den == int: 
                fase_total_den=0
                pass
            fase_total_den += np.degrees(np.angle(termo_den))


        fase_total = fase_total_num-fase_total_den

        diferenca = fase_total - (-180)


        
        if diferenca > 0:
            print(f'Fase total calculada: {round(fase_total,4)}°')
            print(f'Fase para correção: {round(diferenca,4)}°')
        else:
            print(f'Fase total calculada: {round(fase_total,4)}°')
            print(f'Fase para correção: {round(-1*diferenca,4)}°')

    def Condicao_de_magnitude(num_H_s,den_H_s,num_C_s,den_C_s,num_G_s,den_G_s,polo_desejado):

        num_H = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in num_H_s]
        den_H = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in den_H_s]
        num_G = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in num_G_s]
        den_G = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in den_G_s]
        num_C = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in num_C_s]
        den_C = [term.subs(s, polo_desejado) if isinstance(term, sp.Basic) and s in term.free_symbols else term for term in den_C_s]
        num_H= [complex(val) for val in num_H]
        den_H= [complex(val) for val in den_H]
        num_G= [complex(val) for val in num_G]
        den_G= [complex(val) for val in den_G]
        num_C= [complex(val) for val in num_C]
        den_C= [complex(val) for val in den_C]

        numerador_C = 1
        denominador_C = 1
        numerador_H = 1
        denominador_H = 1
        numerador_G = 1
        denominador_G= 1       

        for termo in num_H:
            numerador_C*= termo
        for termo in den_H:
            denominador_C*= termo
        for termo in num_G:
            numerador_G*= termo
        for termo in den_G:
            denominador_G*= termo
        for termo in num_C:
            numerador_H*= termo
        for termo in den_C:
            denominador_H*= termo

        C_s = numerador_C / denominador_C
        G_s = numerador_G / denominador_G
        H_s = numerador_H / denominador_H


        K = 1 / (H_s * G_s * C_s)

        fase = np.degrees(np.angle(K))
        print(f'O valor de K é {round(abs(K), 4)} com fase {round(fase, 4)}°')
    
    def Calculo_T_s(num_H_s,den_H_s,num_C_s,den_C_s,num_G_s,den_G_s):

        numerador_C = 1
        denominador_C = 1
        numerador_H = 1
        denominador_H = 1
        numerador_G = 1
        denominador_G= 1       

        for termo in num_C_s:
            numerador_C*= termo
        for termo in den_C_s:
            denominador_C*= termo
        for termo in num_G_s:
            numerador_G*= termo
        for termo in den_G_s:
            denominador_G*= termo
        for termo in num_H_s:
            numerador_H*= termo
        for termo in den_H_s:
            denominador_H*= termo

        C_s = numerador_C / denominador_C
        G_s = numerador_G / denominador_G
        H_s = numerador_H / denominador_H

        T_s = simplify(C_s * G_s / (1 + C_s * G_s * H_s))

        numerador, denominador = fraction(T_s) 

        denominador = [float(coef) for coef in denominador.as_poly().all_coeffs()]


        raizes = np.roots(denominador)

        print("Função de transferência em malha fechada T(s):", T_s)
        print("\nPolos do sistema em malha fechada:")
        for raiz in raizes:
            print(np.round(raiz, 4))

    def resposta_degrau(num_H_s,den_H_s):
        numerador_H=1
        denominador_H=1

        for termo in num_H_s:
            numerador_H*= termo
        for termo in den_H_s:
            denominador_H*= termo
        
        polinomios_num= sp.expand(numerador_H)
        polinomios_den= sp.expand(denominador_H) 

        num = list(sp.Poly(polinomios_num, s).all_coeffs())
        den = list(sp.Poly(polinomios_den, s).all_coeffs())

        num_ft = [float(val) for val in num]
        den_ft = [float(val) for val in den]

        Func_Transferencia_ft = signal.TransferFunction(num_ft,den_ft)
        t_ft, resposta_ft = signal.step(Func_Transferencia_ft)


        plt.plot(t_ft, resposta_ft)
        plt.title('Resposta ao Degrau - FT')
        plt.xlabel('Tempo [s]')
        plt.ylabel('Amplitude')
        plt.grid(True)

        plt.tight_layout()
        plt.show()

    def Lugar_geometrico_das_raizes(num_H_s,den_H_s):
        numerador_H=1
        denominador_H=1

        for termo in num_H_s:
            numerador_H*= termo
        for termo in den_H_s:
            denominador_H*= termo
        
        polinomios_num= sp.expand(numerador_H)
        polinomios_den= sp.expand(denominador_H) 

        num = list(sp.Poly(polinomios_num, s).all_coeffs())
        den = list(sp.Poly(polinomios_den, s).all_coeffs())

        num_ft = [float(val) for val in num]
        den_ft = [float(val) for val in den]

        Func_Transferencia_ft = tf(num_ft,den_ft)

        # Plot do Lugar Geométrico das Raízes
        plt.figure()
        # Adicione plot=True para evitar problemas de retorno no futuro
        rlocus(Func_Transferencia_ft, plot=True, grid=True)
        plt.title('Lugar geométrico das Raízes - FT')
        plt.xlabel('Real')
        plt.ylabel('Imaginário')
        plt.grid(True)
        plt.show()

    def Diagrama_Bode(num_H_s,den_H_s):
        numerador_H=1
        denominador_H=1

        for termo in num_H_s:
            numerador_H*= termo
        for termo in den_H_s:
            denominador_H*= termo
        
        polinomios_num= sp.expand(numerador_H)
        polinomios_den= sp.expand(denominador_H) 

        num = list(sp.Poly(polinomios_num, s).all_coeffs())
        den = list(sp.Poly(polinomios_den, s).all_coeffs())

        num_ft = [float(val) for val in num]
        den_ft = [float(val) for val in den]

        file_path = r'DADOS.xlsx'
        df = pd.read_excel(file_path)
        df['Frequencia (Rad/s)'] = pd.to_numeric(df['Frequencia (Rad/s)'], errors='coerce')
        df['Magnitude sinal de saida'] = pd.to_numeric(df['Magnitude sinal de saida'], errors='coerce')
        df = df.dropna(subset=['Frequencia (Rad/s)', 'Magnitude sinal de saida'])

        frequencies = df['Frequencia (Rad/s)'].values
        magnitude = df['Magnitude sinal de saida'].values
        fase = df['Defasagem Angular '].values
        magnitude_dB = 20 * np.log10(magnitude)

        fig, axs = plt.subplots(2, 2, figsize=(8, 6))
        ax1, ax2, ax3, ax4 = axs.flatten()

        ax1.plot(frequencies, magnitude_dB, label="Magnitude (dB)")
        ax1.set_xscale('log')
        ax1.set_ylabel("Magnitude (dB)")
        ax1.set_title("Diagrama de magnitude a partir do Excel")
        ax1.grid(True, which="both", ls="--")

        ax2.plot(frequencies, fase, label="Fase (graus)", color='orange')
        ax2.set_xscale('log')
        ax2.set_ylabel("Fase (graus)")
        ax2.set_xlabel("Frequência (Hz)")
        ax2.set_title("Diagrama de fase a partir do Excel")
        ax2.grid(True, which="both", ls="--")

        Func_Transferencia = signal.TransferFunction(num_ft,den_ft) 
        w, mag, fase2 = signal.bode(Func_Transferencia)

        ax3.plot(w, mag, label="Magnitude (dB) FT")
        ax3.set_xscale('log')
        ax3.set_ylabel("Magnitude (dB)")
        ax3.set_title("Diagrama de Bode a partir da H(s)")
        ax3.grid(True, which="both", ls="--")

        ax4.plot(w, fase2, label="Fase (graus)", color='orange')
        ax4.set_xscale('log')
        ax4.set_ylabel("Fase (graus)")
        ax4.set_xlabel("Frequência (Hz)")
        ax4.set_title("Diagrama de fase a partir da H(s)")
        ax4.grid(True, which="both", ls="--")
        plt.tight_layout()
        plt.show()
    
### exemplo de uso #####
##G(s) = 10/s(s-1)(s+1)(s+2)
Polo_desejado = -0.25+ 1j

num_H_s =[10]
den_H_s =[s,s-1,s+1,s+2]

num_C_s=[s + 0.44,s + 0.44]
den_C_s=[s + 5.921]


num_G_s=[1]
den_G_s=[1]
### exemplo de uso #####

controle_classico.Condicao_de_Fase(num_H_s,den_H_s,Polo_desejado)