'''
本代码是对医疗援助后的隔离强度进行敏感性分析
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
plt.rcParams['font.sans-serif'] = ['Arial Unicode MS']
#读取处理好的数据集
data = pd.read_csv('武汉市1.18之后的数据.csv')
prob1 = 0.01671555
prob2 = 0.05222604
prob = (prob1+prob2)/2
gamma = 1/14
d = 0.021
alpha = 1/7 #潜伏期按照7天进行计算
area = 8569.15 #单位：平方千米
people = 1121.2 #单位：万人
people = people * 1e4#单位：人
dist = 0.15 #单位：km
t1 = 1
t2 = 0.3#交通管控后的流动系数
rio = math.sqrt(area) / (math.sqrt(people) - 1)
n1 = math.ceil(t1*dist/(2*rio))**2+4#未防控时的接触人数
n2 = math.ceil(t2*dist/(2*rio))**2+4#交通管控后的接触人数
n3 = 4 #小区管控时的接触人数
k = 2.92815369 #医疗物资系数
l = 0.7#口罩的防护系数
M = 11911#援助前的武汉市总病床数
M1 = 1000+1600+13000 #火神山+雷神山+方舱医院的病床数
f = 0#患者与病床的比例
w_1 = 0.5
w_2 = 0.75
w_3 = 1
w_4 = 1.25
w_5 = 1.5
control_limit = 50 #认为剩余感染者小于等于50即基本控制
def calculate(S,E,I,R,D,w1):
    S[0][0] = people - 1
    E[0][0] = 0
    I[0][0] = 1
    R[0][0] = 0
    D[0][0] = 0
    d1 = 46 #1.22及之前均未进行自我防护(45:12.8~1.22)，d2和d1相等
    d3 = 46+12#2.3及之前均未引入医疗援助（12:1.23～2.3）
    d4 = 46+12+7 #2.10及之前均未引入小区管控(7:2.4~2.10)
    for i in range(0,d1-1):
        f = I[0][i]/M 
        if f <= 1:
            w = 2
        else:
            w = 2+n1*(f-1)
        S[0][i+1] = S[0][i] - (prob * n1)*((E[0][i]*S[0][i])/people) - (prob * w)*((I[0][i]*S[0][i])/people)
        E[0][i+1] = E[0][i] + (prob * n1)*((E[0][i]*S[0][i])/people) + (prob * w)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
        I[0][i+1] = I[0][i] + alpha*E[0][i] - gamma*I[0][i] - d*I[0][i]
        R[0][i+1] = R[0][i] + gamma*I[0][i]
        D[0][i+1] = D[0][i] + d*I[0][i]
    for i in range(d1-1,d3-1):
        f = I[0][i]/M
        if f <= 1:
            w = 2
        else:
            w = 2+n2*(f-1)
        S[0][i+1] = S[0][i] - (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) - (1-l)*(prob * w)*((I[0][i]*S[0][i])/people)
        E[0][i+1] = E[0][i] + (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) + (1-l)*(prob * w)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
        I[0][i+1] = I[0][i] + alpha*E[0][i] - gamma*I[0][i] - d*I[0][i]
        R[0][i+1] = R[0][i] + gamma*I[0][i]
        D[0][i+1] = D[0][i] + d*I[0][i]
    for i in range(d3-1,d4-1):
        f = I[0][i]/(M+M1)
        if f <= 1:
            w = w1
        else:
            w = w1+n2*(f-1)
        S[0][i+1] = S[0][i] - (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) - (1-l)*(prob * w)*((I[0][i]*S[0][i])/people)
        E[0][i+1] = E[0][i] + (1-l)*(prob * n2)*((E[0][i]*S[0][i])/people) + (1-l)*(prob * w)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
        I[0][i+1] = I[0][i] + alpha*E[0][i] - k*gamma*I[0][i] - d*I[0][i]
        R[0][i+1] = R[0][i] + k*gamma*I[0][i]
        D[0][i+1] = D[0][i] + d*I[0][i]
    for i in range(46+12+7-1,365-1):#小区管控
        f = I[0][i]/(M+M1)
        if f <= 1:
            w = w1
        else:
            w = w1+n3*(f-1)
        S[0][i+1] = S[0][i] - (1-l)*(prob * n3)*((E[0][i]*S[0][i])/people) - (1-l)*(prob * w)*((I[0][i]*S[0][i])/people)
        E[0][i+1] = E[0][i] + (1-l)*(prob * n3)*((E[0][i]*S[0][i])/people) + (1-l)*(prob * w)*((I[0][i]*S[0][i])/people) - alpha*E[0][i]
        I[0][i+1] = I[0][i] + alpha*E[0][i] - k*gamma*I[0][i] - d*I[0][i]
        R[0][i+1] = R[0][i] + k*gamma*I[0][i]
        D[0][i+1] = D[0][i] + d*I[0][i]
#w1
S1 = np.zeros((1,365))
E1 = np.zeros((1,365))
I1 = np.zeros((1,365))
R1 = np.zeros((1,365))
D1 = np.zeros((1,365))
#w2
S2 = np.zeros((1,365))
E2 = np.zeros((1,365))
I2 = np.zeros((1,365))
R2 = np.zeros((1,365))
D2 = np.zeros((1,365))
#w3
S3 = np.zeros((1,365))
E3 = np.zeros((1,365))
I3 = np.zeros((1,365))
R3 = np.zeros((1,365))
D3 = np.zeros((1,365))
#w4
S4 = np.zeros((1,365))
E4 = np.zeros((1,365))
I4 = np.zeros((1,365))
R4 = np.zeros((1,365))
D4 = np.zeros((1,365))
#w5
S5 = np.zeros((1,365))
E5 = np.zeros((1,365))
I5 = np.zeros((1,365))
R5 = np.zeros((1,365))
D5 = np.zeros((1,365))
#最终结果
calculate(S1,E1,I1,R1,D1,w_1)
calculate(S2,E2,I2,R2,D2,w_2)
calculate(S3,E3,I3,R3,D3,w_3)
calculate(S4,E4,I4,R4,D4,w_4)
calculate(S5,E5,I5,R5,D5,w_5)
#求对应的指标
max1 = np.max(I1[0,:])
max2 = np.max(I2[0,:])
max3 = np.max(I3[0,:])
max4 = np.max(I4[0,:])
max5 = np.max(I5[0,:])
def control(value,I,E):
    for i in range(len(I[0,:])):
        if I[0][i] == value:
            index = i
            break
    for i in range(len(I[0,:])):
        if I[0][i] <= control_limit and i>index and E[0][i] <= control_limit:
            c = i+1
            break
    return c
c1 = control(max1,I1,E1)
c2 = control(max2,I2,E2)
c3 = control(max3,I3,E3)
c4 = control(max4,I4,E4)
c5 = control(max5,I5,E5)
with open('医疗援助后的隔离效果的敏感性分析.txt','w') as f:
    f.write('援助后的隔离效果为%.2f'%w_1+'最大值为%d，基本控制需要的时间为第%d天\n'%(round(max1),c1))
    f.write('援助后的隔离效果为%.2f'%w_2+'最大值为%d，基本控制需要的时间为第%d天\n'%(round(max2),c2))
    f.write('援助后的隔离效果为%.2f'%w_3+'最大值为%d，基本控制需要的时间为第%d天\n'%(round(max3),c3))
    f.write('援助后的隔离效果为%.2f'%w_4+'最大值为%d，基本控制需要的时间为第%d天\n'%(round(max4),c4))
    f.write('援助后的隔离效果为%.2f'%w_5+'最大值为%d，基本控制需要的时间为第%d天\n'%(round(max5),c5))
f.close()
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(figsize=(8,6))
plt.plot(time[0,:],I1[0,:],color='orange')
plt.plot(time[0,:],I2[0,:],color='hotpink')
plt.plot(time[0,:],I3[0,:],color='red')
plt.plot(time[0,:],I4[0,:],color='skyblue')
plt.plot(time[0,:],I5[0,:],color='lightgreen')
plt.legend(['w1=0.5','w1=0.75','w1=1','w1=1.25','w1=1.5'])
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('医疗援助后的隔离效果对于患者数量的影响.png')