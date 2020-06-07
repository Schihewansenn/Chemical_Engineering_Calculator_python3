# 牛顿迭代法求解函数根
# 有两个模式，求解多项式函数或着任意函数
# 不需要调用第三方库可以直接搞
# 任意函数的求解方法使用说明在最后，多项式函数的可以仿照写

# 差分牛顿二分法
def neSolve(
	D = 1e-12, x0 = 0, f = 0, maxConvergance = 50):
	x = x0
	precision = 1000*D
	for i in range(maxConvergance):
		f1 = testFunc(x)
		f2 = testFunc(x + D)
		k1 = D / (f2 - f1)
		xr = k1 * (f - f1) + x
		# 收敛判断
		if abs(x - xr) < precision :   # 精度符合，收敛，跳出循环
			x = xr
			print('迭代收敛')
			break
		elif abs(f2 - f1) < 5e-15 :
			print('零点畸变')
			break
		elif i < (maxConvergance - 1) : # 若未到最大循环次数且精度不合适，继续迭代
			x = xr		
			print('正在进行第', i + 1,'次迭代，y = ', f1, 'x = ', x)
		else :                         # 意外不收敛
			print('迭代未收敛至要求，结果为: x = ', x)
		
	return x


#获取参数，多项式牛顿二分法迭代求解,返回多项式函数根
def polySolve( 
	maxConvergance = 50,    # 最大迭代次数，默认50
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
	for i in range(maxConvergance) : 
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
		elif i < (maxConvergance - 1) : # 若未到最大循环次数且精度不合适，继续迭代
			x = xr
			print('正在进行第', i + 1,'次迭代，y = ', y, 'x = ', x)
		else :                         # 意外不收敛
			print('迭代未收敛至要求，结果为: x = ', x)

	return x                           # 返回函数值


# 在这里写上你要求的函数
def testFunc(x):
	return 500*(1-1/(1+0.05/12)**x)/(0.05/12)**x

# 在参数种确定差分精度、迭代初值、目标函数值、最大收敛次数(顺序相同)
a = neSolve(D = 1e-13, x0 = 1.5, f = 1000000, maxConvergance = 100)
print('函数根 = ', a)