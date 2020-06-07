'----------------------------cutline---------------------------------'
# stage-by-stage equilibrium calculation model
'----------------------------cutline---------------------------------'
import time
import os
import matplotlib.pyplot as plt

# 函数区

# 汽液相平衡方程（气相组成，相对挥发度）：返回液相组成
def phaseEquilibrium(y, alpha):
	a1 = 1 / alpha
	a2 = 1 / y
	a3 = 1 / (a1 + a2 - 1)
	x = a1 * a3
	return x

# 精馏段操作线方程（本级液相组成，回流比，塔顶组成）：返回下一级气相组成
def rectifyingOperating(x, R, xD):
	a1 = R / (R+1)
	a2 = xD / (R+1)
	y = a1*x + a2
	return y

# 精馏段操作线方程（本级液相组成，降液采出比，塔顶组成）：返回下一级气相组成
# 降液采出比：Rt = L' / W
def strippingOperating(x, Rt, xW):
	a1 = Rt / (Rt - 1)
	a2 = xW / (Rt - 1) 
	y = a1*x - a2
	return y

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

# 降液采出比函数，输入回流比与塔顶馏出比，返回塔釜降液采出比
def rTCAL(R, DF):
	a1 = 1 + R*DF
	a2 = 1 - DF
	a3 = a1 / a2
	return a3

# 塔顶采出比计算函数，输入进料组成、规定的塔顶组成和塔釜组成，返回塔顶采出比
def dfCAL(xF, xD, xW):
	a1 = xF - xW
	a2 = xD - xW
	a3 = a1 / a2
	return a3

# 最小回流比计算函数，输入塔顶组成、进料组成组成和相对挥发度，返回塔顶采出比
def minR(xD, xF, alpha):
	ye = xF*alpha/(1+(alpha-1)*xF)
	a1 = xD - xF
	a2 = ye - xF
	a3 = (a1 / a2) - 1
	return a3

def cutline(words = 'cutline', lenth = 40, symble = '_'):
    a1 = len(words)
    a2 = (lenth - 2*a1) / 2
    a3 = round(a2)
    print("\n", end = '')

    for i in range(a3) : 
        print(symble, end = '')

    print(words, end = '')

    for j in range(a3):
        print(symble, end = '')

    print("\n \n", end = '')

'----------------------------cutline---------------------------------'

cutline(words = '模块说明')
print('本模块为板式塔的精馏塔设计型计算模块，默认泡点进料与泡点回流，使用全凝器，使用逐板计算算法。')
print('输入组成信息：塔顶、塔釜、进料；回流比（或倍数）')
print('输出组成信息：各层塔板气/液相数据；最小回流比；理论塔板数；塔顶采出比')
print('未特殊说明各单位均按mol计算，其他采取IS制')
cutline(words = '请按照要求输入信息')

# 获取参数
alpha = eval(input("请输入相对挥发度: "))
xF = eval(input("请输入进料组成(molFrac) xF: "))
xD = eval(input("请输入塔顶组成(molFrac) xD: "))
xW = eval(input("请输入塔釜组成(molFrac) xW: "))
Rmin = minR(xD, xF, alpha)
print('最小回流比 Rmin = ', Rmin)
R = eval(input("请输入回流比(> 0)或最小回流比的倍数(< -1)(mole base): "))

if R > 0 :
	R = R
elif R < 0 :
	R = -R*Rmin
else:
	print('输入回流比数值错误')

Eff = eval(input("请输入全塔效率 (0 < eff <= 1): "))

cutline(words = '开始初始化')

# 初始化
Flag = prChack(Eff)
Accuracy = 6
xD = round(xD, Accuracy)
x = 1
y = xD
N = 0
M = 0
NE = 0
xf = 0
yf = 0
DF = dfCAL(xF, xD, xW)
Rt = rTCAL(R, DF)
Convergance = 1000
Rectx = []
Strix = []
Recty = [y]
Striy = []
Totalx = []
Totaly = []

'----------------------------cutline---------------------------------'
cutline(words = '开始计算')

if Flag == False :
	# 条件循环逐板计算，核心函数

	# 精馏段计算循环
	for i in range(Convergance):               # 最大收敛次数                        
		x = phaseEquilibrium(y, alpha)         # 汽-液相平衡方程
		if prChack(x):                         # 精馏段液相组成数据检查
			break                              # 若有错误则跳出循环
		
		if x > xF :                            # 判断是否是精馏段
			x = round(x, Accuracy)             # 保留指定位数
			Rectx.append(x)                    # 将精馏段各板液相数据组成列表
			y = rectifyingOperating(x, R, xD)  # 精馏段操作线方程
			if prChack(y):                     # 精馏段液相组成数据检查
				break                          # 若有错误则跳出循环
			y = round(y, Accuracy)
			Recty.append(y)                    # 将精馏段各板汽相数据组成列表
			N = N + 1                          # 精馏段塔板数加一
			print("第", N,"块精馏段塔板计算完毕 \n")
		
		else :                                 # 将最后一块塔板液相组成存入数组
			Recty.pop(-1)
			print("精馏段计算结束 \n ")   
			break                              
	if N == Convergance :
		print('精馏段计算不收敛 \n ')

	# 二次初始化
	y = Recty[-1]
	X = Rectx[-1]

	# 进料板组成计算
	yf = round(strippingOperating(x, Rt, xW), Accuracy)
	xf = round(phaseEquilibrium(yf, alpha), Accuracy)

	# 三次初始化
	x = xf
	y = yf

	# 提馏段计算循环
	for i in range(Convergance):               # 最大收敛次数                        
		
		if x > xW :                            # 判断是否是精馏段
			y = strippingOperating(x, Rt, xW)  # 提馏段操作线方程
			if prChack(y):                     # 提馏段液相组成数据检查
				break                          # 若有错误则跳出循环
			y = round(y, Accuracy)             # 将数值保留至指定精度
			Striy.append(y)                    # 将提馏段各板汽相数据组成列表
			x = phaseEquilibrium(y, alpha)     # 汽-液相平衡方程
			if prChack(x):                     # 提馏段液相组成数据检查
				break                          # 若有错误则跳出循环
			x = round(x, Accuracy)             # 保留指定位数
			Strix.append(x)                    # 将提馏段各板液相数据组成列表
			M = M + 1                          # 提馏段塔板数加一
			print("第", M,"块提馏段塔板计算完毕 \n")
		
		else :
			print("提馏段计算结束 \n")          # 计算终止
			break                              # 终止循环

	if M == Convergance :
		print('提馏段计算不收敛 \n')

else : 
	print('输入组成有误，请检查，程序终止请重试')


# 数据分布调整

Totalx = Rectx + [xf] + Strix
Totaly = Recty + [yf] + Striy


#实际塔板数计算
NE = N + M + 1
Nt = NE - 1

NR = Nt / Eff         # 计算实际塔板数
NR = round(NR)        # 取整
N2 = N / Eff          # 计算实际进料板位置
N2 = round(N2)        # 取整
'----------------------------cutline---------------------------------'
#计算结果输出

cutline(words = "计算概况")

print("理论板数NE(含再沸器与分凝器) = ", NE, '\n')
print("精馏段塔板数NR = ", N-1, '\n')
print("提馏段塔板数NS = ", M, '\n')
print("进料板位置N = ", N, '\n')
print("实际塔板数Nt = ", NR, '\n')
print("实际塔板进料位置Nrt = ", N2, '\n')
print('塔顶馏出比(mole) = ', DF, '\n')
print('降液采出比(mole) = ', Rt, '\n')

cutline(words = "组成信息")

print("精馏段液相组成 = ", Rectx, '\n')
print("精馏段气相组成 = ", Recty, '\n') 
print("进料板液相组成 = ", xf, '\n')
print("进料板气相组成 = ", yf, '\n')
print("提馏段液相组成 = ", Strix, '\n')
print("提馏段气相组成 = ", Striy, '\n')
print("全塔液相组成 = ", Totalx, '\n')
print("全塔气相组成 = ", Totaly, '\n')

plt.title('x(y)-NE')
plt.xlim(1, len(Totalx)+1)
plt.ylim(0, 1)
plt.plot(range(1, len(Totalx)+1), Totalx, '-', color = 'b',label = 'x')  
plt.plot(range(1, len(Totaly)+1), Totaly, '-', color = 'r', label = 'y')
plt.show()

cutline(words = "输出完成，计算成功")
