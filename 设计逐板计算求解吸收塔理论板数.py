'----------------------------cutline---------------------------------'
# stage-by-stage equilibrium calculation model
# 轻组分为吸收质
# 重组分为吸收剂
# 设计计算主程序
# 校核计算主程序
'----------------------------cutline---------------------------------'
import time
import os
import matplotlib.pyplot as plt
import numpy as np 


# 函数区

# 汽液相平衡方程（气相组成，相对挥发度）：返回液相组成
def phaseEquilibrium(y, m):
    a1 = 1 / m
    a2 = 1 / y
    a3 = 1 / (a1 + a2 - 1)
    x = a1 * a3
    return x

# 操作线方程（本级液相组成，上级液相组成，本级气相组成）：返回下级气相组成
def Operating(x1, x0, y1, LV):
    a1 = x1 - x0
    y2 = a1*LV + y1
    return y2

# 气相分率转化函数（气相分压，吸收总压）：返回气相分率
def moleFracy(p, pt):
    return p / pt

# 数据校核函数，判断分数是否超出定义域，如果出现错误返回Ture，否则返回False
def prChack(Frac):   
    if Frac < 0 :
        print("工艺参数错误，汽/液相组成小于0")
        return True
    elif Frac > 1 :
        print("工艺参数错误，汽/液相组成大于1")
        return True
    else :
        return False

def cutline(words = 'cutline', lenth = 40, symble = '_'):
    a1 = len(words)
    a2 = (lenth - a1) / 2
    a3 = round(a2)
    print("\n", end = '')

    for i in range(a3) : 
        print(symble, end = '')

    print(words, end = '')

    for j in range(a3):
        print(symble, end = '')

    print("\n \n", end = '')


'----------------------------cutline---------------------------------'

cutline(words = '请按照要求输入信息')

# 获取参数
m = eval(input("亨利常数 m: "))
yW = eval(input("请输入进入气体组成(molFrac) yw: "))
y1 = eval(input("请输入输出气体组成(molFrac)上限 y1: "))
x0 = eval(input("请输入新鲜吸收剂组成(molFrac) x0: "))
LV = eval(input("请输入液气比(mole base): "))
Eff = eval(input("请输入全塔效率 (0 < eff <= 1): "))

cutline(words = '开始初始化')

# 初始化
Accuracy = 6
x = x0
y = y1
N = 0
NE = 0
Convergance = 100
xlst = [x]
ylst = [y]
xls = []
yls = []

'----------------------------cutline---------------------------------'
cutline(words = '开始计算')

# 条件循环逐板计算，核心函数

# 设计型计算循环
for i in range(Convergance):                               # 最大收敛次数                        
    x = phaseEquilibrium(y, m)                             # 汽-液相平衡方程
    if prChack(x):                                         # 精馏段液相组成数据检查
        break                                              # 若有错误则跳出循环
    
    if y < yW :                                            # 判断进气组成上限是否小于实际进气组成
        x = round(x, Accuracy)                             # 保留指定位数
        xlst.append(x)                                    # 将精馏段各板液相数据组成列表
        y = Operating(xlst[-1], xlst[-2], ylst[-1], LV)    # 精馏段操作线方程
        if prChack(y):                                     # 精馏段液相组成数据检查
            break                                          # 若有错误则跳出循环
        y = round(y, Accuracy)
        ylst.append(y)                                    # 将精馏段各板汽相数据组成列表
        N = N + 1                                          # 精馏段塔板数加一
        print("第", N,"块塔板计算完毕 \n")
    
    else :                                                 # 将最后一块塔板液相组成存入数组
        x0 = xlst[0]
        yW = ylst[-1]
        ylst.pop(-1)
        xlst.pop(0)
        N = N + 1
        print("计算结束 \n ")   
        break                              
if N == Convergance :
    print('计算不收敛 \n ')


NR = N / Eff         # 计算实际塔板数         


'----------------------------cutline---------------------------------'
#计算结果输出

cutline(words = "计算概况")

print("理论板数N = ", N, '\n')
print("实际塔板数Nt = ", NR, '\n')

cutline(words = "组成信息")

print("液相组成 = ", xlst, '\n')
print("气相组成 = ", ylst, '\n') 
print("吸收液组成 = ", xlst[-1], '\n')
print("进气组成 = ", yW, '\n')

plt.title('塔内气/液相组成-理论塔板数分布')
plt.xlim(1, len(xlst)+1)
plt.ylim(0, yW*1.1)
plt.plot(range(1, len(xlst)+1), xlst, '-.', color = 'b',label = 'x')  
plt.plot(range(1, len(ylst)+1), ylst, '-*', color = 'r', label = 'y')
plt.show()

cutline(words = "输出完成，计算成功")
