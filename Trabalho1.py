import pandas as pd
import numpy as np
import math
import cmath
import matplotlib.pyplot as plt
import tkinter as tk
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg


root = tk.Tk()
root.title('Gráfico Interativo')

H_value = tk.StringVar()
P_value = tk.StringVar()
Scc_value = tk.StringVar()
Tcc_value = tk.StringVar()


H_label = tk.Label(root, textvariable=H_value)
H_label.pack()

P_label = tk.Label(root, textvariable=P_value)
P_label.pack()

Scc_label = tk.Label(root, textvariable=Scc_value)
Scc_label.pack()

Tcc_label = tk.Label(root, textvariable=Tcc_value)
Tcc_label.pack()

canvas = None


def plotar_grafico(P, H):

    global canvas, H_value, P_value, Scc_value, Tcc_value

    if canvas is not None:
        canvas.get_tk_widget().destroy()

    Q = complex(0, 0.074)
    S = P + Q
    U = 1
    F = 60
    delta_t = 0.01

    # 2º Passo: Obtenção da corrente I
    I = S.conjugate()/U

    # 3º Passo: Obtenção da reatância equivalente
    X_antes = complex(0, 0.65)
    X_durante = complex(0, 1.8)
    X_depois = complex(0, 0.8)

    # 4º Passo: Obtenção da fem E
    E = U + I*X_antes
    Em, EO = cmath.polar(E)
    EO = EO * 180/math.pi

    # 5º Passo: Obtenção de Pmáx
    Pm_antes = Em*U/abs(X_antes)
    Pm_durante = Em*U/abs(X_durante)
    Pm_depois = Em*U/abs(X_depois)

    # 6º Passo: Obtenção de So
    So = math.asin(P/Pm_antes)

    # 7º Passo: Obtenção de Smax
    Sf = math.asin(P/Pm_depois)
    Smax = math.pi - Sf

    # 8º Passo: Obtenção de Scc
    Scc = math.acos((P*So - P*Smax - Pm_depois*math.cos(Smax) +
                    Pm_durante*math.cos(So))/(Pm_durante - Pm_depois))
    Scc_ang = Scc*180/math.pi

    # 9º Passo: Obtenção do Tcc
    K = 180*F*(delta_t**2)/H

    Pe = [0, 0, 0, 0, 0, 0]
    Pa = [0, 0, 0, 0, 0, 0]
    KPa = [0, 0, 0, 0, 0, 0]
    deltaS = [0, 0, 0, 0, 0, 0]
    Si = [So, 0, 0, 0, 0, 0]
    lista = [2, 3, 4]
    Pe[0] = Pm_antes*math.sin(Si[0])
    Pe[1] = Pm_durante*math.sin(Si[0])
    Pa[0] = P - Pe[1]
    Pa[1] = (0 + Pa[0])/2
    KPa[1] = Pa[1]*K
    deltaS[1] = deltaS[0] + KPa[1]
    Si[1] = So

    for i in lista:
        Pe[i] = Pm_durante*math.sin(Si[i-1]+(deltaS[i-1]*math.pi/180))
        Pa[i] = P - Pe[i]
        KPa[i] = Pa[i]*K
        deltaS[i] = deltaS[i-1] + KPa[i]
        Si[i] = Si[i-1] + (deltaS[i]*math.pi/180)

    matriz_angulos = np.array(
        [[Si[2]*180/math.pi], [Si[3]*180/math.pi], [Si[4]*180/math.pi]])
    matriz_quadrada = np.array(
        [[1, 0.01, 0.0001], [1, 0.02, 0.0004], [1, 0.03, 0.0009]])
    matriz_inversa = np.linalg.inv(matriz_quadrada)
    matriz_a = np.dot(matriz_inversa, matriz_angulos)
    c, b, a = matriz_a[0], matriz_a[1], matriz_a[2]
    c = c - Scc*180/math.pi
    Tcc_m = (-b + ((b**2)-4*a*c)**0.5)/(2*a)
    Tcc = Tcc_m[0]

    # Atualizar os valores das variáveis
    H_value.set(f'H: {H:.2f}')
    P_value.set(f'Pm: {P:.2f}')
    Scc_value.set(f'Scc: {(Scc*180/math.pi):.2f}')
    Tcc_value.set(f'Tcc: {Tcc:.3f}')

    Se = np.linspace(0, np.pi, 180)
    Pe_antes = Pm_antes*np.sin(Se)
    Pe_durante = Pm_durante*np.sin(Se)
    Pe_depois = Pm_depois*np.sin(Se)
    Pem = np.full_like(Se, P)

    fig, ax = plt.subplots()
    ax.plot(Se, Pe_antes, label='Antes da falta')
    ax.plot(Se, Pe_durante, label='Durante a falta')
    ax.plot(Se, Pe_depois, label='Depois da falta')
    ax.plot(Se, Pem, label='Pm')

    ax.fill_between(Se, Pem, Pe_durante, where=(
        (Se >= So) & (Se <= Scc)), color='gray', alpha=0.5)
    ax.fill_between(Se, Pe_depois, Pem, where=(
        (Se >= Scc) & (Se <= Smax)), color='gray', alpha=0.5)

    ax.set_xlabel('Ângulo')
    ax.set_ylabel('Potência')
    ax.set_title('Gráfico do ângulo de chaveamento crítico')
    ax.legend()
    canvas = FigureCanvasTkAgg(fig, master=root)
    canvas.draw()
    canvas.get_tk_widget().pack()


def atualizar_grafico():
    try:
        P = float(pm_entry.get())
        H = float(H_entry.get())
        plotar_grafico(P, H)
    except ValueError:
        print("Erro: Certifique-se de que os valores de P e H estão corretos.")


pm_label = tk.Label(root, text='Pm:')
pm_label.pack()
pm_entry = tk.Entry(root)
pm_entry.pack()

H_label = tk.Label(root, text='H:')
H_label.pack()
H_entry = tk.Entry(root)
H_entry.pack()


atualizar_btn = tk.Button(
    root, text='Atualizar Gráfico', command=atualizar_grafico)
atualizar_btn.pack()


def ao_fechar():
    root.destroy()
    root.quit()

root.protocol("WM_DELETE_WINDOW", ao_fechar)


root.mainloop()
