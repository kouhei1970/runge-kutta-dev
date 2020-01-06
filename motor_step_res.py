import matplotlib.pyplot as plt
import numpy as np

'''
def rk4(func, t, h, y, *x)
ルンゲ・クッタ法を一回分計算する関数
    引数リスト
    func:導関数
    t：現在時刻を表す変数
    h：刻み幅
    y：出力変数（求めたい値）
    *x:引数の数が可変する事に対応する、その他の必要変数
※この関数では時刻は更新されないため、これとは別に時間更新をする必要があります。
'''
def rk4(func, t, h, y, *x):
    k1=h*func(t, y, *x)
    k2=h*func(t+0.5*h, y+0.5*k1, *x)
    k3=h*func(t+0.5*h, y+0.5*k2, *x) 
    k4=h*func(t+h, y+k3, *x)
    y=y+(k1 + 2*k2 + 2*k3 + k4)/6
    return y

'''
導関数の書き方
def func(t, y, *state):
    func:自分で好きな関数名をつけられます
    t:時刻変数(変数の文字はtで無くても良い) 
    y:出力変数(変数の文字はyで無くても良い)
    *state:その他の必要変数(引数の数は可変可能))
#関数サンプル
def vdot(t, y, *state):
    s1=state[0]
    s2=state[1]
    return t+y+s1+s2
    
'''

def idot(t, i, omega, e):
    R=1.07
    K=1.98e-3
    L=17e-6
    return (e - R*i -K*omega)/L 

def omegadot(t, omega, i , T_L):
    K=1.98e-3
    D=1.226e-7
    J=0.59e-7
    return (K*i - D*omega - T_L)/J

'''
SecondOrderDelaySysStepRes(zeta, omega, Kg, t)
２次振動系のステップ応答の解析解を計算する
パラメータ（引数）
zeta:減衰係数
omega:自然各周波数
Kg:ゲイン
t:時間配列（計算時刻のnumpy配列）
'''

def SecondOrderDelaySysStepRes(zeta, omega, Kg, t):

    if zeta > 1:

        y=Kg * (1 - np.exp(-zeta*omega*t) * 
                 (    (zeta + np.sqrt(zeta**2 - 1) )*np.exp( omega*np.sqrt(zeta**2 - 1)*t) 
                    - (zeta - np.sqrt(zeta**2 - 1) )*np.exp(-omega*np.sqrt(zeta**2 - 1)*t)
                 ) /2/np.sqrt(zeta**2 - 1)
               )
    elif zeta < 1:
        phi=np.arctan2(np.sqrt(1-zeta**2), zeta)
        y=Kg * ( 1 - np.exp(-zeta*omega*t) * np.sin(omega*np.sqrt(1 - zeta**2)*t + phi)
                / np.sqrt(1 - zeta**2)
               )
    else:
        y=Kg*(1 - np.exp(-omega*t)*(omega*t  + 1 ))

    return y
'''

２次振動系の係数からゲイン、減衰係数、自然各周波数を求める

         D
G(s)=---------
     As^2+Bs+C
から
              K omega^2
G(s)=-----------------------------
     s^2 + 2 zeta omega s +omega^2

K(ゲイン),zeta(減衰係数),omega（自然各周波数）を計算する。
     
'''
def find2ndOrderDelaySysParameter(A,B,C,D):
    K=D/C
    zeta=B/2/np.sqrt(A*C)
    omega=np.sqrt(C/A)
    return zeta, omega, K


#ここからメイン
def main():
    #初期化
    t=0.0
    omega=0.0
    i=0.0
    e=3.0
    h=1.0e-6
    Omega=[]
    I=[]
    T=[]
    endtime=0.02
    Num=int(endtime/h)+1

    #求解ループ
    for n in range(Num):
        Omega.append(omega)
        I.append(i)
        T.append(t)
        omegaold=omega
        iold=i
        i=rk4(idot, t, h, iold, omegaold, e)
        omega=rk4(omegadot, t, h, omegaold, iold, 0)
        t=t+h
    #解析解計算
    #モータのパラメータ
    R=1.07
    K=1.98e-3
    L=17e-6   
    D=1.226e-7
    J=0.59e-7
    
    A=J*L
    B=D*L + J*R
    C=D*R + K**2
    #print(A,B,C) #２次振動系ステップ応答解析解関数がオーバーフローするのでデバグ用に追加
    
    #２次振動系のパラメータ計算
    zeta, omega_n, Kg = find2ndOrderDelaySysParameter(A,B,C,K*3)
    #print(zeta, omega_n, Kg) #上記関数がオーバーフローしたのでデバグ用に追加
    #計算時刻配列作成
    t=np.linspace(0, endtime, 1000)
    #解析解計算
    y=SecondOrderDelaySysStepRes(zeta, omega_n, Kg, t)
    
    plt.plot(T, Omega, label='Numerical')
    plt.plot(t, y, label='Analytical')
    plt.xlabel('Time(s)')
    plt.ylabel('omega(rad/s)')
    plt.grid()
    plt.legend()
    plt.show()

if __name__=='__main__':
    main()
