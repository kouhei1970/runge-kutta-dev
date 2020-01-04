'''
作りはじめ日時　2019/12/31
ルンゲクッタ法のを使うシミュレーションプログラムのための開発
バネマスダンパ系のシムレーションで検証
どう書けば各種シミュレーションに使いやすいのかを
検討していきたい。
'''

#%matplotlib inline
import numpy as np
import matplotlib.pyplot as plt

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

#バネマスダンパ系速度の導関数
def vdot(t, v, *state):
    D=3
    K=10
    M=1.0
    x=state[0]
    u=state[1]
    return (-D*v -K*x +u)/M

#バネマスダンパ系位置の導関数
def xdot(t, x, *state):
    v=state[0]
    return v

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


#ここからメイン
def main():
    #ルンゲ・クッタ法による数値計算 
    #初期化
    t=0.0
    v=0.0
    x=0.0
    u=1.0
    h=0.001
    V=[]
    X=[]
    T=[]

    #求解ループ
    for n in range(5000):
        V.append(v)
        X.append(x)
        T.append(t)
        vold=v
        xold=x
        v=rk4(vdot, t, h, vold, xold, u)
        x=rk4(xdot, t, h, xold, vold)
        t=t+h
    
    #解析解に基づく計算
    M=1.0
    D=3.0
    Ks=10.0
        
    #Kg=2.0
    #zeta=1.1
    #omega=2*np.pi
    Kg=1/Ks
    zeta=D/2/np.sqrt(Ks*M)
    omega=np.sqrt(Ks/M)

    print('Kg={:.3f} zeta={:.3f} Omega_n={:.3f}'.format(Kg, zeta, omega))
    
    t=np.linspace(0, 5, 1000)

    y = SecondOrderDelaySysStepRes(zeta, omega, Kg, t)

    plt.figure(figsize=(5, 3))
    plt.plot(T, X)
    plt.plot(t, y)
    plt.grid()
    plt.show()
    

if __name__=='__main__':
    main()
