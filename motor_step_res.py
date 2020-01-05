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



#ここからメイン
def main():
    #初期化
    t=0.0
    omega=0.0
    i=0.0
    e=3.0
    h=1.0e-7
    Omega=[]
    I=[]
    T=[]
    endtime=0.15
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

    plt.plot(T, Omega)
    plt.xlabel('Time(s)')
    plt.ylabel('omega(rad/s)')
    plt.grid()
    plt.show()

if __name__=='__main__':
    main()

