'------------ThermodynamicPack------------'
# 热力学计算包
# 数值参数输入函数

# 多项式求解函数（4阶，牛顿迭代法）

# RK方程a、b参数计算函数
# 求压缩因子RK方程（RKV）
# 求压力RK方程（RKp）
# 求温度RK方程（RKT）

# 安托因方程

# wilson 相平衡模型，公制单位，温度是K，压力 帕（不完善，不建议使用）
# 对偶参数计算函数
# 活度系数计算函数
# 理想气体逸度计算函数
# 不可压缩液体逸度计算函数
# 求解气相组成函数
# 求解液相组成函数
# 求解总压计算函数
'-----------------------------------------'
import numpy as np

# 分割线函数，让结果井然有序
def cutline(words = 'cutline', lenth = 40, symble = '-'):
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


#获取参数（输入提示词），返回数值
def cDia(words) :
	return eval(input(words))

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
	b = RKb(pc, Tc, R)	                         # 计算参数 b
	dia1 = R*T /(V-b)                            
	dia2 = a / (pow(T, 0.5)*V*(V+b))             
	p = dia1 - dia2                              # 计算压强 p
	return p                                     # 返回压强 p


# 求解温度的RK方程（临界压力，临界温度，温度，压力，摩尔气体常数），返回温度T（K）
def RKT(pc, Tc, p, V, R = 8.314) :
	a = RKa (pc, Tc, R)                          # 计算参数 a
	b = RKb(pc, Tc, R)	                         # 计算参数 b
	a0 = -pow(((V-b)/(V+b))*(a/V), 2)            # 温度迭代方程常数项
	a1 = pow(p*(V-b), 2)                         # 温度迭代方程一次项系数
	a2 = -2*p*R*(V-b)                            # 温度迭代方程二次项系数
	a3 = pow(R, 2)                               # 温度迭代方程三次项系数
	
	for i in range (10) : # 调用求解函数
		T = polySolve(x0 = 50*i, a0 = a0, a1 = a1, a2 = a2, a3 = a3)
		if 150 < T : # 判断计算结果是否符合实验规律，是否更换初值计算
			break
	
	return T                                      # 返回温度 T

# 对偶参数计算方程（第一组分液相摩尔体积、第二组分液相摩尔体积、相互作用参数、温度、气体常数），返回对偶参数
def compareDiameter(V1L, V2L, G, T, R = 8.314):
	a1 = V2L / V1L
	a2 = G / (R*T)
	A12 = a1 * np.exp(-a2) 
	return A12

# 活度系数计算方程（第一组分液相摩尔体积、第二组分液相摩尔体积、12相互作用参数、21相互作用参数、第一组分组成、温度）
def gamma(V1L, V2L, G12, G21, x1, T) :
	x2 = 1 - x1
	a1 = compareDiameter(V1L, V2L, G12, T)
	a2 = compareDiameter(V2L, V1L, G21, T)
	a3 = -np.log(x1 + a1*x2)
	a4 = (a1 / (x1 + a1*x2) - (a2 / (x2 + a2*x1)))
	a5 = a3 + x2*a4 
	a6 = -np.log(x2 + a2*x1)
	a7 = a6 - x1*a4
	gamma1 = np.exp(a5) 
	gamma2 = np.exp(a7)
	return gamma1, gamma2

# Antoin方程
def steamPressure(A, B, C, T) :
	a1 = A - B / ( T + C )
	ps = pow(10, a1) * 1000
	return ps

# 不可压缩液体逸度
def incompressibleFL(gamma, x, ps) :
	return gamma * x * ps

# 理想气体逸度
def idealFV(y, p) :
	return y * p	

# 相平衡计算气相分率
def gasFrac(x, gamma, p, ps) :
	return (x * gamma * ps) / p

# 相平衡计算液相分率
def liqFrac(y, gamma, p, ps) :
	return (y * p) / (gamma * ps)

# 相平衡计算压力
def wilsonP(y, gamma, x, ps) :
	return (gamma * x * ps) / y

'------------------------------------------------------------------------------'
def main() :
	cutline('欢迎使用热力学计算器')
	print('1.牛顿迭代法求解多项式根')
	print('2.RK方程计算压力')
	print('3.RK方程计算压缩因子')
	print('4.RK方程计算温度')
	print('5.Antoin方程计算饱和蒸汽压')
	print('6.Wilson活度系数模型计算活度系数 \n ')
	par = input('请输入计算编号: ')
	cutline('开始计算，所有单位无特殊说明均按照SI制')
	if par == '1' :
		print('正在计算多项式的根')
		return polySolve( 
		maxCycloTimes = 50,                       # 最大迭代次数，默认50
		precision = 1e-13,                        # 绝对精度
		x0 = cDia('迭代初值 = '),                  # 迭代初值，默认0
		f = 0,                                    # 目标函数值，默认0
		a0 = cDia('常数项系数 = '),                # 多项式常数项，默认0
		a1 = cDia('一次项系数 = '),                # 多项式一次项系数，默认0
		a2 = cDia('二次项系数 = '),                # 多项式二次项系数，默认0
		a3 = cDia('三次项系数 = '),                # 多项式三次项系数，默认0
		a4 = cDia('四次项系数 = ')                 # 多项式四次项系数，默认0
		) 
	elif par == '2' :
		print('正在使用RK方程计算压力')
		return RKP(
		pc = cDia('临界压力 = '), 
		Tc = cDia('临界温度 = '), 
		V = cDia('摩尔体积 = '), 
		T = cDia('温度 = '), 
		R = 8.314
		)
	elif par == '3' :
		print('正在使用RK方程计算压缩因子')
		return RKV(
		pc = cDia('临界压力 = '), 
		Tc = cDia('临界温度 = '), 
		p = cDia('压力 = '), 
		T = cDia('温度 = '), 
		R = 8.314
		)
	elif par == '4' :
		print('正在使用RK方程计算温度')
		return RKT(
		pc = cDia('临界压力 = '), 
		Tc = cDia('临界温度 = '), 
		p = cDia('压力 = '), 
		V = cDia('摩尔体积 = '), 
		R = 8.314
		)
	elif par == '5' :
		print('正在使用Antoin方程计算饱和蒸汽压/常用 kPa')
		return steamPressure(
		A = cDia('A = '), 
		B = cDia('B = '), 
		C = cDia('C = '), 
		T = cDia('温度 = ')
		)
	elif par == '6' :
		print('正在使用Wilson方程计算1组分活度系数')
		return gamma(
		V1L = cDia('1组分饱和液体体积 = '), 
		V2L = cDia('2组分液体饱和体积 = '), 
		G12 = cDia('12组分相互作用能量参数 = '), 
		G21 = cDia('21组分相互作用能量参数 = '), 
		x1 = cDia('组分1摩尔分数 = '), 
		T = cDia('温度 = ')
		)
	elif par == '7' :
		print('正在建设中')
	elif par == '8' :
		print('正在建设中')
	elif par == '9' :
		print('正在建设中')
	else :
		print('错误')
		
print(main())
while True :
	if input('继续使用热力学计算器请输入(Y):') == 'Y':
		print(main())
	else :
		break