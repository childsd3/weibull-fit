#!/usr/bin/env python3
import math
import statistics
import glob
import matplotlib.pyplot as plt
import numpy as np

def calculateXXSum (xList):
	xxSum = 0
	for item in xList:
		xxSum += item*item
	return xxSum

def calculateXYSum (xList, yList):
	if len(xList) != len(yList):
		return 'lists dont match!'
	xySum = 0
	for i in range(len(xList)):
		xySum += xList[i]*yList[i]
	return xySum

def calculateListSum (inputList):
	listSum = 0
	for item in inputList:
		listSum += item
	return listSum

def calculateXPrimeList (xList):
	outputList = list()
	for item in xList:
		outputList.append(math.log(item))
	return outputList

def calculateYPrimeList (yList):
	outputList = list()
	for item in yList:
		outputList.append(math.log(math.log(1/(1-item))))
	return outputList

def shiftXValues (xList,alpha):
	outputList = list()
	for item in xList:
		outputList.append(item - alpha)
	return outputList

def calculateXPrimes (xList):
	outputList = list()
	for item in xList:
		outputList.append(math.log(item))
	return outputList

def calculateCDF (listLength):
	outputList = list()
	for item in range(1,listLength+1):
		outputList.append(item/(listLength+1))
	return outputList

def calculateWeibullMean (alpha,beta,m):
	return math.pow(m,1/beta)*math.gamma(1+1/beta)+alpha

def calculateWeibullStandardDeviation (alpha,beta,m):
	return math.sqrt(math.pow(m,2/beta)*(math.gamma(1+2/beta)-math.pow(math.gamma(1+1/beta),2)))

def generateWeibullList (alpha,beta,m):
	return alpha

def calculateRegression (filename,alpha):
	values = []
	with open(filename, "r") as f:
		for line in f:
			if float(line):
				values.append(float(line))
	values.sort()
	if values[0] <= alpha:
		print('Alpha needs to be less than the smallest data point:',values[0])
		return
	n = len(values)
	shiftedXValues = shiftXValues(values,alpha)
	xPrime = calculateXPrimeList(shiftedXValues)
	cdfValues = calculateCDF(n)
	yPrime = calculateYPrimeList(cdfValues)
	xxSum = calculateXXSum(xPrime)
	xySum = calculateXYSum(xPrime,yPrime)
	xSum = calculateListSum(xPrime)
	ySum = calculateListSum(yPrime)
	averageXPrime = statistics.mean(xPrime)
	averageYPrime = statistics.mean(yPrime)
	beta = (n*xySum-xSum*ySum)/(n*xxSum-math.pow(xSum,2))
	m = math.exp(beta*averageXPrime-averageYPrime)
	weibullMean = calculateWeibullMean(alpha,beta,m)
	weibullStandardDeviation = calculateWeibullStandardDeviation(alpha,beta,m)
	print('mean: '+str(weibullMean))
	return {
		'alpha': alpha,
		'beta': beta,
		'm': m,
		'weibullMean': weibullMean,
		'weibullStandardDeviation': weibullStandardDeviation,
		'values': values
	}

for item in glob.glob('*.dat'):
	regression = calculateRegression(item,0)
	print('******************************************************')
	print('                  filename: '+item)
	print('                     alpha: '+str(regression['alpha']))
	print('                      beta: '+str(regression['beta']))
	print('                         m: '+str(regression['m']))
	print('              Weibull Mean: '+str(regression['weibullMean']))
	print('Weibull Standard Deviation: '+str(regression['weibullStandardDeviation']))
	print('******************************************************')
#	hist, bin_edges = np.histogram(regression['values'], density=True, bins='auto')
	def oneMinus (value):
		return 1 - value
	def weibullPlot (x,alpha,beta,m):
		return math.exp(math.pow(x-alpha,beta)/(-m))
	def negExpPlot (x,average):
		return math.exp(-x/average)
	oneMinus = np.vectorize(oneMinus)
	weibullPlot = np.vectorize(weibullPlot)
	negExpPlot = np.vectorize(negExpPlot)
	regressionLength = len(regression['values'])+1
	CDF = np.arange(1/regressionLength,1,1/regressionLength)
	oneMinusValues = oneMinus(CDF)
	weibullValues = weibullPlot(regression['values'], regression['alpha'], regression['beta'], regression['m'])
	negExpValues = negExpPlot(regression['values'],np.average(regression['values']))
	plt.plot(regression['values'],oneMinusValues, label='Measured')
	plt.plot(regression['values'],weibullValues, label='Fit')
	plt.plot(regression['values'],negExpValues, label='Negative Exponential')
	plt.legend()
	plt.show()