'''
作りはじめ日時　2019/12/31

ルンゲクッタ法のを使うシミュレーションプログラムのための開発
バネマスダンパ系のシムレーションで検証

どうかけば各種シミュレーションに使いやすいのかを
検討していきたい。

'''

#%matplotlib inline
import matplotlib.pyplot as plt

def rk4(func, t, h, y, *x):
    k1=h*func(t, y, *x)
    k2=h*func(t+0.5*h, y+0.5*k1, *x)
    k3=h*func(t+0.5*h, y+0.5*k2, *x) 
    k4=h*func(t+h, y+k3, *x)
    y=y+(k1 + 2*k2 + 2*k3 + k4)/6
    #print(y)    
    return y

'''

導関数の書き方
def func(t, y, *state):

func(時刻変数, 出力変数, その他の必要変数(引数の数は可変可能))

#関数サンプル
def func(t, y, *state):
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

t=0.0
v=0.0
x=0.0
u=1.0
h=0.01
V=[]
X=[]
T=[]
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


