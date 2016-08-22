#!/usr/bin/env python3
import math
import statistics
import argparse

parser = argparse.ArgumentParser(description='CNI Weibull Fitting Program')
parser.add_argument('INPUT')
parser.add_argument('ALPHA')
parser.add_argument('-p','--plot', action='store_true', help='Use matplotlib to display results')
parser.add_argument('-c','--csv', action='store_true', help='Export results as CSV')

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
#	print('              alpha:  '+str(alpha))
#	print('               beta:  '+str(beta))
#	print('                  m:  '+str(m))
#	print('               Mean:  '+str(weibullMean))
#	print(' Standard Deviation:  '+str(weibullStandardDeviation))
	return (alpha,beta,m,weibullMean,weibullStandardDeviation)

calculateRegression(parser.parse_args().INPUT,float(parser.parse_args().ALPHA))

if parser.parse_args().csv:
	print(calculateRegression(parser.parse_args().INPUT,float(parser.parse_args().ALPHA)))

#if parser.parse_args().plot:
#	import matplotlib
#	print('plot!')