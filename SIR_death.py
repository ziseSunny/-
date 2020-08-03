'''
本代码是在basic版本的基础上考虑了死亡率，对模型进行适当修正
'''
import numpy as np
import pandas as pd
import math
import datetime
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
#读取处理好的数据集
data = pd.read_csv('武汉市1.18之后的数据.csv')
class estimate:
    def __init__(self,C,timeIndex,n,gamma,d,x0,steps):
        self.timeRange = np.array([i for i in range(timeIndex[0], timeIndex[1]+1)])
        self.n = n
        self.gamma = gamma
        self.d = d
        self.x0 = x0
        self.steps = steps 
        self.C = C
    def residual(self,prob):#对应的残差函数
        res = np.array(self.C*np.exp((prob * self.n - self.d - self.gamma) * self.timeRange) - data.loc[self.timeRange - self.steps,'累计确诊'])
        return (res**2).sum() / self.timeRange.size
    def optimize(self):
        self.solution = minimize(self.residual, self.x0, method='nelder-mead', options={'xtol':1e-8, 'disp':True})
        print('Infection Probability: ',self.solution.x)
        return self.solution.x
    def getBasicReproduct(self):
        self.basicRep = self.n * self.solution.x[0] / (self.gamma+self.d)
        print('Basic Reproduction Number: ',self.basicRep)
        return self.basicRep
startTime = datetime.datetime.strptime('2019-12-08', "%Y-%m-%d")
timeBox = [datetime.datetime.strptime('2020-01-18', "%Y-%m-%d"), datetime.datetime.strptime('2020-01-22', "%Y-%m-%d")]
timeIndex = [(t - startTime).days for t in timeBox]
#假设人群在武汉市均匀分布，求取每个人的平均接触人数
area = 8569.15 #单位：平方千米
people = 1121.2 #单位：万人
#由于病毒传播还是依靠于较近距离的接触，只计算每个人的接触人数，设每个人的每日集中活动距离为0.15km
avg_dist = math.sqrt(area) / (math.sqrt(people*1e4) - 1)
dist = 0.075
#假定圆形范围内的人均满足要求
num = int(dist / avg_dist)
#根据给定的距离生成所需的矩阵，进行判断
flags = np.zeros((2*num+1,2*num+1))
for i in range(0,2*num+1):
    for j in range(0,2*num+1):
        if (i-num)**2+(j-num)**2 <= num*num:
            flags[i][j] = 1
#统计数量
n= 0
for i in range(len(flags)):
    for j in range(len(flags[i])):
        if flags[i][j]:
            n = n+1
gamma = 1/14
d = 0.021
x0 = 0.02
C = 1
steps = 41#从12.8~1.18共41天
estimation = estimate(C, timeIndex, n, gamma, d, x0, steps)
infect_prob = estimation.optimize()
beta = n*infect_prob
basicReproduction = estimation.getBasicReproduct()
R0=beta / (gamma+d)
#提取第二阶段的传播系数（1.24～2.3）
startTime_2 = datetime.datetime.strptime('2020-01-24', "%Y-%m-%d")
timeBox_2 = [datetime.datetime.strptime('2020-01-24', "%Y-%m-%d"), datetime.datetime.strptime('2020-02-04', "%Y-%m-%d")]
timeIndex_2 = [(t - startTime_2).days for t in timeBox_2]
n_2 = 4
x0_2=0.02
C_2 = 0#1月24日武汉市对应的累计确诊病例
d_2 = d
gamma_2 = gamma
steps_2 = 0
for i in range(len(data.index)):
    if data.loc[i,'报道时间'] == '1月24日':
        C_2 = data.loc[i,'累计确诊']
        break
estimation_2 = estimate(C_2, timeIndex_2, n_2, gamma_2, d_2, x0_2, steps_2)
infect_prob2 = estimation_2.optimize()
beta_2 = n_2*infect_prob2
basicReproduction_2 = estimation_2.getBasicReproduct()
R02 = beta_2 / (gamma_2+d_2)
#对新冠病毒传播疫情的两阶段模拟（仍然将观察时间定为1年）
S = np.zeros((1,365))
I = np.zeros((1,365))
R = np.zeros((1,365))
D = np.zeros((1,365))
S[0][0] = people*1e4 - 1
I[0][0] = 1
R[0][0] = 0
D[0][0] = 0
for i in range(0,46-1):
    S[0][i+1] = S[0][i] - (beta * I[0][i] * S[0][i])/(people*1e4)
    I[0][i+1] = I[0][i] + (beta * I[0][i] * S[0][i])/(people*1e4) - gamma * I[0][i] - d * I[0][i]
    R[0][i+1] = R[0][i] + gamma*I[0][i]
    D[0][i+1] = D[0][i] + d*I[0][i]
for i in range(46-1,365-1):
    S[0][i+1] = S[0][i] - (beta_2 * I[0][i] * S[0][i])/(people*1e4)
    I[0][i+1] = I[0][i] + (beta_2 * I[0][i] * S[0][i])/(people*1e4) - gamma_2 * I[0][i] - d_2 * I[0][i]
    R[0][i+1] = R[0][i] + gamma_2*I[0][i]
    D[0][i+1] = D[0][i] + d_2*I[0][i]
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(figsize=(8, 6)) 
plt.plot(time[0,:],I[0,:],color='red')
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('改进后的SIR模型对应的患者人数趋势图.png')
max_I = np.max(I[0,:])
index_max = 0
control_limit = 50
for i in range(len(I[0,:])):
    if I[0][i] == max_I:
        index_max = i
print('最大患者数:%d，对应的时间：%d'%(round(max_I),index_max))
for i in range(len(I[0,:])):
    if I[0][i] <= control_limit and i > index_max:
        index = i
        break
print('基本控制住的时间：%d'%index)