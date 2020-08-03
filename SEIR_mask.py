'''
本代码是在SEIR模型基础上探讨口罩因素的引入,假设居民在封城后才意识到戴口罩
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
#读取处理好的数据集
data = pd.read_csv('武汉市1.18之后的数据.csv')
prob1 = 0.01671555
prob2 = 0.05222604
prob = (prob1+prob2)/2
gamma = 1/14
d = 0.021
alpha = 1/7 #潜伏期按照7天进行计算
n1 = 13 #未防控时的接触人数
n2 = 4 #防控时的接触人数
k = 2.92815369 #医疗物资系数
l = 0.7
S = np.zeros((1,365))
E = np.zeros((1,365))
I = np.zeros((1,365))
R = np.zeros((1,365))
D = np.zeros((1,365))
people = 1121.2 * 1e4 #单位：人
S[0][0] = people - 1
E[0][0] = 0
I[0][0] = 1
R[0][0] = 0
D[0][0] = 0
for i in range(0,46-1):
    S[0][i+1] = S[0][i] - (prob * n1)*((E[0][i]*S[0][i])/people) - (prob * n1)*((I[0][i]*S[0][i])/people)
    E[0][i+1] = E[0][i] + (prob * n1)*((E[0][i]*S[0][i])/people) + (prob * n1)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
    I[0][i+1] = I[0][i] + alpha*E[0][i] - gamma*I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + gamma*I[0][i]
    D[0][i+1] = D[0][i] + d*I[0][i]
for i in range(46-1,46+12-1):
    S[0][i+1] = S[0][i] - (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) - (1-l)*(prob * n2)*((I[0][i]*S[0][i])/people)
    E[0][i+1] = E[0][i] + (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) + (1-l)*(prob * n2)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
    I[0][i+1] = I[0][i] + alpha*E[0][i] - gamma*I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + gamma*I[0][i]
    D[0][i+1] = D[0][i] + d*I[0][i]
for i in range(46+12-1,365-1):
    S[0][i+1] = S[0][i] - (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) - (1-l)*(prob * n2)*((I[0][i]*S[0][i])/people)
    E[0][i+1] = E[0][i] + (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) + (1-l)*(prob * n2)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
    I[0][i+1] = I[0][i] + alpha*E[0][i] - k*gamma*I[0][i] - d*I[0][i]
    R[0][i+1] = R[0][i] + k*gamma*I[0][i]
    D[0][i+1] = D[0][i] + d*I[0][i]
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(figsize=(8, 6)) 
plt.plot(time[0,:],I[0,:],color='red')
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('引入口罩后的SEIR模型对应的患者人数趋势图.png')
#求解疫情基本得到控制的时间和患病总人数
control_limit = 50
max_I = np.max(I[0,:])
index_max = 0
index = 0
for i in range(len(I[0,:])):
    if I[0][i] == max_I:
        index_max = i
print('最大患者数:%d，对应的时间：%d'%(round(max_I),index_max))
for i in range(0,365):
    if I[0][i] <= control_limit and i > index_max and E[0][i] <= control_limit:
        index = i+1
        break
print('基本控制住的时间：%d'%index)