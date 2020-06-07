'------------ThermodynamicPack------------'
# 热力学计算包
# 数值参数输入函数
# 多项式求解函数（4阶）
# RK方程a、b参数计算函数
# 求压缩因子RK方程（RKV）
# 求压力RK方程（RKp）
# 求温度RK方程（RKT）
'-----------------------------------------'
import numpy as np
import time

#获取参数（输入提示词），返回数值
def catchDiameters(words) :
    return eval(input(words, '(单位: K, KPa, L)'))

#获取参数，牛顿二分法迭代求解,返回多项式函数根
def polySolve( 
    maxCycloTimes = 50,    # 最大迭代次数，默认50
    precision = 1e-13,     # 绝对精度
    x0 = 0,                # 迭代初值，默认0
    f = 0,                 # 目标函数值，默认0
    a0 = 0,                # 多项式常数项，默认0
    a1 = 0,                # 多项式一次项系数，默认0
    a2 = 0,                # 多项式二次项系数，默认0
    a3 = 0,                # 多项式三次项系数，默认0
    a4 = 0                 # 多项式四次项系数，默认0
    ) :

    x = x0
    print('开始方程求解')   # 调用函数提示
    # 牛顿二分法核心算法
    for i in range(maxCycloTimes) : 
        # 计算函数值
        y = a0 + a1*x + a2*pow(x, 2) + a3*pow(x, 3) + a4*pow(x, 4)
        # 计算导数值
        Dy = a1 + 2*a2*x + 3*a3*pow(x, 2) + 4*a4*pow(x, 3)
        # 计算收敛系数（导数倒数）
        k = 1 / Dy
        # 计算切线与目标函数值交点横坐标
        xr = k * (f - y) + x
        # 收敛判断
        if abs(x - xr) < precision :   # 精度符合，收敛，跳出循环
            x = xr
            print('迭代收敛')
            break
        elif i < (maxCycloTimes - 1) : # 若未到最大循环次数且精度不合适，继续迭代
            x = xr
            print('正在进行第', i + 1,'次迭代，y = ', y, 'x = ', x)
        else :                         # 意外不收敛
            print('迭代未收敛至要求，结果为: x = ', x)

    return x                           # 返回函数值
    


# RK方程参数 a（临界压力，临界温度，气体常数）
def RKa(pc, Tc, R):
    a = 0.42748 * pow(R, 2) * pow(Tc, 2.5) / pc
    return a

# RK方程参数 b（临界压力，临界温度，气体常数）
def RKb(pc, Tc, R):
    b = 0.08664 * R * Tc / pc
    return b


# 求解摩尔体积的RK方程（临界压力，临界温度，温度，压力，摩尔气体常数），返回压缩因子Z
def RKV(pc, Tc, p, T, R = 8.314) :               
    a = RKa (pc, Tc, R)                          # 计算参数 a
    b = RKb(pc, Tc,R)                            # 计算参数 b
    A = a*p / (pow(R, 2) * pow(T, 2.5))          # 计算参数 A
    B = b*p / (R*T)                              # 计算参数 B
    a0 = -A*B                                    # 压缩因子迭代方程常数项
    a1 = A - B - pow(B, 2)                       # 压缩因子迭代方程一次项系数
    a2 = 1                                       # 压缩因子迭代方程二次项系数
    a3 = 1                                       # 压缩因子迭代方程三次项系数

    for i in range (10) :      # 调用求解函数
        Z = polySolve(x0 = 0.1*i, a0 = a0, a1 = a1, a2 = a2, a3 = a3)
        if 0 < Z :             # 判断计算结果是否符合实验规律，是否更换初值计算
            break

    return Z                                     # 返回压缩因子 Z


# 求解压强的RK方程（临界压力，临界温度，温度，压力，摩尔气体常数），返回压强p（Pa）
def RKP(pc, Tc, V, T, R = 8.314) :               
    a = RKa (pc, Tc, R)                          # 计算参数 a
    b = RKb(pc, Tc, R)                           # 计算参数 b
    dia1 = R*T /(V-b)                            
    dia2 = a / (pow(T, 0.5)*V*(V+b))             
    p = dia1 - dia2                              # 计算压强 p
    return p                                     # 返回压强 p


# 求解温度的RK方程（临界压力，临界温度，温度，压力，摩尔气体常数），返回温度T（K）
def RKT(pc, Tc, p, V, R = 8.314) :
    a = RKa (pc, Tc, R)                          # 计算参数 a
    b = RKb(pc, Tc, R)                           # 计算参数 b
    a0 = -pow(((V-b)/(V+b))*(a/V), 2)            # 温度迭代方程常数项
    a1 = pow(p*(V-b), 2)                         # 温度迭代方程一次项系数
    a2 = -2*p*R*(V-b)                            # 温度迭代方程二次项系数
    a3 = pow(R, 2)                               # 温度迭代方程三次项系数
    
    for i in range (10) : # 调用求解函数
        T = polySolve(x0 = 50*i, a0 = a0, a1 = a1, a2 = a2, a3 = a3)
        if 150 < T : # 判断计算结果是否符合实验规律，是否更换初值计算
            break
    
    return T                                      # 返回温度 T


t1 = time.time()
print(RKV(3.65e6, 408.1, 3.704e5, 300))
print(RKP(3.39e6, 126.2, 0.04636e-3, 273.15))
print(RKT(3.648e6, 408.1, 3.704e5, 6.081e-3))
print(RKV(3.65e6, 408.1, 3.704e5, 300))
print(RKP(3.39e6, 126.2, 0.04636e-3, 273.15))
print(RKT(3.648e6, 408.1, 3.704e5, 6.081e-3))
print(RKV(3.65e6, 408.1, 3.704e5, 300))
print(RKP(3.39e6, 126.2, 0.04636e-3, 273.15))
print(RKT(3.648e6, 408.1, 3.704e5, 6.081e-3))
print(RKV(3.65e6, 408.1, 3.704e5, 300))
print(RKP(3.39e6, 126.2, 0.04636e-3, 273.15))
print(RKT(3.648e6, 408.1, 3.704e5, 6.081e-3))
print(RKV(3.65e6, 408.1, 3.704e5, 300))
print(RKP(3.39e6, 126.2, 0.04636e-3, 273.15))
print(RKT(3.648e6, 408.1, 3.704e5, 6.081e-3))
t2 = time.time()
print(t2-t1,'s')



