'--------------------------------------------------'
# wilson 相平衡模型，公制单位，温度是K，压力 帕
'--------------------------------------------------'

import numpy as np

def compareDiameter(V1L, V2L, G, T, R = 8.314):
	a1 = V2L / V1L
	a2 = G / (R*T)
	A12 = a1 * np.exp(-a2) 
	return A12


def gamma(V1L, V2L, G12, G21, x1, T):
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

def steamPressure(A, B, C, T):
	a1 = A - B / ( T + C )
	ps = pow(10, a1) * 1000
	return ps


def incompressibleFL(gamma, x, ps):
	return gamma * x * ps

def idealFV(y, p):
	return y * p	

def gasFrac(x, gamma, p, ps):
	return (x * gamma * ps) / p

def liqFrac(y, gamma, p, ps):
	return (y * p) / (gamma * ps)

def wilsonP(y, gamma, x, ps):
	return (gamma * x * ps) / y

def wilsonT(y, gamma, x, ps):
	return (gamma * x * ps) / y