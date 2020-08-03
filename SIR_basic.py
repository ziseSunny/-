'''
本代码是利用最基本的SIR模型对武汉市的新冠肺炎疫情传播进行模拟
'''
import numpy as np
import pandas as pd
import math
import datetime
from scipy.optimize import minimize
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
dataset = pd.read_csv('Updates_NC.csv') #读取整个数据集
#从中提取出武汉市的数据
data_total = dataset[dataset['城市']=='武汉市']
#观察得到的数据，进行累加操作
Infection = data_total.groupby('报道时间')['新增确诊'].sum()#这样处理考虑到一天可能有多个数据
Recovered = data_total.groupby('报道时间')['新增出院'].sum()
Dead = data_total.groupby('报道时间')['新增死亡'].sum()
data_process = {'报道时间':Infection.index, '新增确诊':Infection.values, '新增出院':Recovered.values, '新增死亡':Dead.values}
#将data_process转换为pd的DataFrame格式
data_process = pd.DataFrame(data_process, index=[i for i in range(Infection.shape[0])])
#进行累加
Infect_Total = [data_process.loc[0:i,'新增确诊'].sum() for i in range(data_process.shape[0])]
Recover_Total = [data_process.loc[0:i,'新增出院'].sum() for i in range(data_process.shape[0])]
Dead_Total = [data_process.loc[0:i,'新增死亡'].sum() for i in range(data_process.shape[0])]
#将累加的数据添加进表格中
data_process = data_process.join(pd.DataFrame([Infect_Total,Recover_Total,Dead_Total],index=['累计确诊','累计出院','累计死亡']).T)#.T表示表格转置
#删除1.18之前的数据
data = data_process[data_process['报道时间']>='1月18日']
data.index = [i for i in range(data.shape[0])]#重新编辑表格的编号
#对传染率进行预测，由于1.23采取封城措施，只考虑1.22及之前的数据，由于12.8~1.17的数据存在缺失，选取1.18~1.22的数据进行模拟
class estimate:
    def __init__(self,timeIndex,n,gamma,x0):
        self.timeRange = np.array([i for i in range(timeIndex[0], timeIndex[1]+1)])
        self.n = n
        self.gamma = gamma
        self.x0 = x0
        self.steps = 41 #从12.8~1.18共41天(此时的12.8为第0天，1。18为第41天)
    def residual(self,prob):#对应的残差函数
        res = np.array(np.exp((prob * self.n - self.gamma) * self.timeRange) - data.loc[self.timeRange - self.steps,'累计确诊'])
        return (res**2).sum() / self.timeRange.size
    def optimize(self):
        self.solution = minimize(self.residual, self.x0, method='nelder-mead', options={'xtol':1e-8, 'disp':True})
        print('Infection Probability: ',self.solution.x)
        return self.solution.x
    def getBasicReproduct(self):
        self.basicRep = self.n * self.solution.x[0] / (self.gamma)
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
dist = 0.075#人的活动半径
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
x0 = 0.02
estimation = estimate(timeIndex, n, gamma, x0)
infect_prob = estimation.optimize()
beta = n*infect_prob
basicReproduction = estimation.getBasicReproduct()
R0=beta / gamma
#模拟数据
I = np.zeros((1,365))
R = np.zeros((1,365))
S = np.zeros((1,365))
I[0][0] = 1
R[0][0] = 0
S[0][0] = people*1e4 - 1
for i in range(0,365-1):
    S[0][i+1] = S[0][i] - (beta * I[0][i] * S[0][i])/(people*1e4)
    I[0][i+1] = I[0][i] + (beta * I[0][i] * S[0][i])/(people*1e4) - gamma * I[0][i]
    R[0][i+1] = R[0][i] + gamma*I[0][i]
#绘制对应的曲线
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(figsize=(8, 6)) 
plt.plot(time[0,:],I[0,:],color='red')
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('最基本的SIR模型对应的患者人数趋势图.png')
control_limit = 50
max_I = np.max(I[0,:])
index_max = 0
for i in range(len(I[0,:])):
    if I[0][i] == max_I:
        index_max = i
print('最大患者数:%d，对应的时间：%d'%(round(max_I),index_max))
for i in range(len(I[0,:])):
    if I[0][i] <= control_limit and i > index_max:
        index = i
        break
print('基本控制住的时间：%d'%index)
data.to_csv('武汉市1.18之后的数据.csv')