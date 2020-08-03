'''
该代码进一步考虑医疗设备的增加带来的改变
'''
import numpy as np
import pandas as pd
import datetime
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
#读取处理好的数据集
data = pd.read_csv('武汉市1.18之后的数据.csv')
class estimate:
    def __init__(self,C,timeIndex,gamma,d,x0,beta):
        self.timeRange = np.array([i for i in range(timeIndex[0], timeIndex[1]+1)])
        self.gamma = gamma
        self.d = d
        self.x0 = x0
        self.C = C
        self.beta = beta
    def residual(self,k):#对应的残差函数
        res = np.array(self.C*np.exp((self.beta - self.d - k*self.gamma) * self.timeRange) - data.loc[self.timeRange,'累计确诊'])
        return (res**2).sum() / self.timeRange.size
    def optimize(self):
        self.solution = minimize(self.residual, self.x0, method='nelder-mead', options={'xtol':1e-8, 'disp':True})
        print('k: ',self.solution.x)
        return self.solution.x
    def getBasicReproduct(self):
        self.basicRep = self.beta / (self.solution.x[0]*self.gamma + self.d)
        print('Basic Reproduction Number: ',self.basicRep)
        return self.basicRep
beta_1 = 0.21730209
beta_2 = 0.20890417
C_1 = 1
C_2 = 0#1月24日武汉市对应的累计确诊病例
for i in range(len(data.index)):
    if data.loc[i,'报道时间'] == '1月24日':
        C_2 = data.loc[i,'累计确诊']
        break
C_3 = 0#1月24日武汉市对应的累计确诊病例
for i in range(len(data.index)):
    if data.loc[i,'报道时间'] == '2月4日':
        C_3 = data.loc[i,'累计确诊']
        break
gamma = 1/14
d = 0.021
x0 = 1
startTime = datetime.datetime.strptime('2020-02-04', "%Y-%m-%d")
timeBox = [datetime.datetime.strptime('2020-02-04', "%Y-%m-%d"), datetime.datetime.strptime('2020-03-10', "%Y-%m-%d")]
timeIndex = [(t - startTime).days for t in timeBox]
estimation = estimate(C_3, timeIndex, gamma, d, x0, beta_2)
k = estimation.optimize()
basicReproduct=estimation.getBasicReproduct()
R0 = beta_2 / (k*gamma+d)
#对数据进行模拟
people = 1121.2 #单位：万人
S = np.zeros((1,365))
I = np.zeros((1,365))
R = np.zeros((1,365))
D = np.zeros((1,365))
S[0][0] = people*1e4 - 1
I[0][0] = 1
R[0][0] = 0
D[0][0] = 0
for i in range(0,46-1):
    S[0][i+1] = S[0][i] - (beta_1 * S[0][i] * I[0][i]) / (people*1e4)
    I[0][i+1] = I[0][i] + (beta_1 * S[0][i] * I[0][i]) / (people*1e4) - gamma * I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + gamma * I[0][i]
    D[0][i+1] = D[0][i] + d * I[0][i]
for i in range(46-1,46+12-1):
    S[0][i+1] = S[0][i] - (beta_2 * S[0][i] * I[0][i]) / (people*1e4)
    I[0][i+1] = I[0][i] + (beta_2 * S[0][i] * I[0][i]) / (people*1e4) - gamma * I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + gamma * I[0][i]
    D[0][i+1] = D[0][i] + d * I[0][i]
for i in range(46+12-1,365-1):
    S[0][i+1] = S[0][i] - (beta_2 * S[0][i] * I[0][i]) / (people*1e4)
    I[0][i+1] = I[0][i] + (beta_2 * S[0][i] * I[0][i]) / (people*1e4) - k*gamma * I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + k*gamma * I[0][i]
    D[0][i+1] = D[0][i] + d * I[0][i]
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(figsize=(8, 6)) 
plt.plot(time[0,:],I[0,:],color='red')
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('引入医疗后的SIR模型对应的患者人数趋势图.png')
control_limit = 50
max_I = np.max(I[0,:])
index_max = 0
index = 0
for i in range(len(I[0,:])):
    if I[0][i] == max_I:
        index_max = i
print('最大患者数:%d，对应的时间：%d'%(round(max_I),index_max))
for i in range(len(I[0,:])):
    if I[0][i] <= control_limit and i > index_max:
        index = i
        break
print('基本控制住的时间：%d'%index)