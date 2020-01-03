'''
作りはじめ日時　2019/12/31

ルンゲクッタ法のを使うシミュレーションプログラムのための開発
バネマスダンパ系のシムレーションで検証

どう書けば各種シミュレーションに使いやすいのかを
検討していきたい。

'''

#%matplotlib inline
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
def xdot(t, x):
    return v

#ここからメイン

#初期化
t=0.0
v=0.0
x=0.0
u=1.0
h=0.01
V=[]
X=[]
T=[]

#求解ループ
for n in range(500):
    V.append(v)
    X.append(x)
    T.append(t)
    v=rk4(vdot, t, h, v, x, u)
    x=rk4(xdot, t, h, x)
    t=t+h

plt.plot(T, X)
plt.grid()
plt.show()


