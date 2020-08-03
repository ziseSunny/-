'''
本代码是探讨不断加入因素对于感染人数变化的影响（因素：无措施、社交距离控制、医疗设施、自我防控、隔离）
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
w = 0#患者感染系数
control_limit = 50 #认为剩余感染者小于等于50即基本控制
#无举措模拟
S0 = np.zeros((1,365))
E0 = np.zeros((1,365))
I0 = np.zeros((1,365))
R0 = np.zeros((1,365))
D0 = np.zeros((1,365))
S0[0][0] = people - 1
E0[0][0] = 0
I0[0][0] = 1
R0[0][0] = 0
D0[0][0] = 0
for i in range(0,365-1):
    S0[0][i+1] = S0[0][i] - (prob * n1)*((E0[0][i]*S0[0][i])/people) - (prob * n1)*((I0[0][i]*S0[0][i])/people)
    E0[0][i+1] = E0[0][i] + (prob * n1)*((E0[0][i]*S0[0][i])/people) + (prob * n1)*((I0[0][i]*S0[0][i])/people) - alpha*E0[0][i]
    I0[0][i+1] = I0[0][i] + alpha*E0[0][i] - gamma*I0[0][i] - d*I0[0][i]
    R0[0][i+1] = R0[0][i] + gamma*I0[0][i]
    D0[0][i+1] = D0[0][i] + d*I0[0][i]
max_I0 = np.max(I0[0,:])
with open('各因素叠加的探讨结果.txt','w') as f:
    f.write('无任何措施的结果：\n')
    f.write('最大值：%d\n'%round(max_I0))
    for i in range(len(I0[0,:])):
        if I0[0][i] == max_I0:
            index_max = i
    for i in range(0,365):
        if I0[0][i] <= control_limit and i > index_max and E0[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n\n'%index)
f.close()
#只进行社交距离控制
S1 = np.zeros((1,365))
E1 = np.zeros((1,365))
I1 = np.zeros((1,365))
R1 = np.zeros((1,365))
D1 = np.zeros((1,365))
#此时认为n为n2,直接认为是封城带来的影响
S1[0][0] = people - 1
E1[0][0] = 0
I1[0][0] = 1
R1[0][0] = 0
D1[0][0] = 0
for i in range(0,46-1):
    S1[0][i+1] = S1[0][i] - (prob * n1)*((E1[0][i]*S1[0][i])/people) - (prob * n1)*((I1[0][i]*S1[0][i])/people)
    E1[0][i+1] = E1[0][i] + (prob * n1)*((E1[0][i]*S1[0][i])/people) + (prob * n1)*((I1[0][i]*S1[0][i])/people) - alpha*E1[0][i]
    I1[0][i+1] = I1[0][i] + alpha*E1[0][i] - gamma*I1[0][i] - d*I1[0][i]
    R1[0][i+1] = R1[0][i] + gamma*I1[0][i]
    D1[0][i+1] = D1[0][i] + d*I1[0][i]
for i in range(46-1,365-1):
    S1[0][i+1] = S1[0][i] - (prob * n2)*((E1[0][i]*S1[0][i])/people) - (prob * n2)*((I1[0][i]*S1[0][i])/people)
    E1[0][i+1] = E1[0][i] + (prob * n2)*((E1[0][i]*S1[0][i])/people) + (prob * n2)*((I1[0][i]*S1[0][i])/people) - alpha*E1[0][i]
    I1[0][i+1] = I1[0][i] + alpha*E1[0][i] - gamma*I1[0][i] - d*I1[0][i]
    R1[0][i+1] = R1[0][i] + gamma*I1[0][i]
    D1[0][i+1] = D1[0][i] + d*I1[0][i]
max_I1 = np.max(I1[0,:])
with open('各因素叠加的探讨结果.txt','a+') as f:
    f.write('调节社交距离的结果：\n')
    f.write('最大值：%d\n'%round(max_I1))
    for i in range(len(I1[0,:])):
        if I1[0][i] == max_I1:
            index_max = i
    for i in range(0,365):
        if I1[0][i] <= control_limit and i > index_max and E1[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n'%index)
    r = (round(max_I0) - round(max_I1))/round(max_I0)*100
    f.write('该措施的提升能力(百分比)：%.2f\n\n'%r)
f.close()
#社交距离+医疗举措
S2 = np.zeros((1,365))
E2 = np.zeros((1,365))
I2 = np.zeros((1,365))
R2 = np.zeros((1,365))
D2 = np.zeros((1,365))
#此时认为n为n2,直接认为是封城带来的影响
S2[0][0] = people - 1
E2[0][0] = 0
I2[0][0] = 1
R2[0][0] = 0
D2[0][0] = 0
for i in range(0,46-1):
    S2[0][i+1] = S2[0][i] - (prob * n1)*((E2[0][i]*S2[0][i])/people) - (prob * n1)*((I2[0][i]*S2[0][i])/people)
    E2[0][i+1] = E2[0][i] + (prob * n1)*((E2[0][i]*S2[0][i])/people) + (prob * n1)*((I2[0][i]*S2[0][i])/people) - alpha*E2[0][i]
    I2[0][i+1] = I2[0][i] + alpha*E2[0][i] - gamma*I2[0][i] - d*I2[0][i]
    R2[0][i+1] = R2[0][i] + gamma*I2[0][i]
    D2[0][i+1] = D2[0][i] + d*I2[0][i]
for i in range(46-1,46+12-1):
    S2[0][i+1] = S2[0][i] - (prob * n2)*((E2[0][i]*S2[0][i])/people) - (prob * n2)*((I2[0][i]*S2[0][i])/people)
    E2[0][i+1] = E2[0][i] + (prob * n2)*((E2[0][i]*S2[0][i])/people) + (prob * n2)*((I2[0][i]*S2[0][i])/people) - alpha*E2[0][i]
    I2[0][i+1] = I2[0][i] + alpha*E2[0][i] - gamma*I2[0][i] - d*I2[0][i]
    R2[0][i+1] = R2[0][i] + gamma*I2[0][i]
    D2[0][i+1] = D2[0][i] + d*I2[0][i]
for i in range(46+12-1,365-1):
    S2[0][i+1] = S2[0][i] - (prob * n2)*((E2[0][i]*S2[0][i])/people) - (prob * n2)*((I2[0][i]*S2[0][i])/people)
    E2[0][i+1] = E2[0][i] + (prob * n2)*((E2[0][i]*S2[0][i])/people) + (prob * n2)*((I2[0][i]*S2[0][i])/people) - alpha*E2[0][i]
    I2[0][i+1] = I2[0][i] + alpha*E2[0][i] - k*gamma*I2[0][i] - d*I2[0][i]
    R2[0][i+1] = R2[0][i] + k*gamma*I2[0][i]
    D2[0][i+1] = D2[0][i] + d*I2[0][i]
max_I2 = np.max(I2[0,:])
with open('各因素叠加的探讨结果.txt','a+') as f:
    f.write('调节社交距离并引入医疗因素的结果：\n')
    f.write('最大值：%d\n'%round(max_I2))
    for i in range(len(I2[0,:])):
        if I2[0][i] == max_I2:
            index_max = i
    for i in range(0,365):
        if I2[0][i] <= control_limit and i > index_max and E2[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n'%index)
    r = (round(max_I1) - round(max_I2))/round(max_I1)*100
    f.write('该措施的提升能力(百分比)：%.2f\n\n'%r)
f.close()
#社交距离+医疗举措+自我防护
S3 = np.zeros((1,365))
E3 = np.zeros((1,365))
I3 = np.zeros((1,365))
R3 = np.zeros((1,365))
D3 = np.zeros((1,365))
#此时认为n为n2,直接认为是封城带来的影响
S3[0][0] = people - 1
E3[0][0] = 0
I3[0][0] = 1
R3[0][0] = 0
D3[0][0] = 0
for i in range(0,46-1):
    S3[0][i+1] = S3[0][i] - (prob * n1)*((E3[0][i]*S3[0][i])/people) - (prob * n1)*((I3[0][i]*S3[0][i])/people)
    E3[0][i+1] = E3[0][i] + (prob * n1)*((E3[0][i]*S3[0][i])/people) + (prob * n1)*((I3[0][i]*S3[0][i])/people) - alpha*E3[0][i]
    I3[0][i+1] = I3[0][i] + alpha*E3[0][i] - gamma*I3[0][i] - d*I3[0][i]
    R3[0][i+1] = R3[0][i] + gamma*I3[0][i]
    D3[0][i+1] = D3[0][i] + d*I3[0][i]
for i in range(46-1,46+12-1):
    S3[0][i+1] = S3[0][i] - (1-l)*(prob * n2)*((E3[0][i]*S3[0][i])/people) - (1-l)*(prob * n2)*((I3[0][i]*S3[0][i])/people)
    E3[0][i+1] = E3[0][i] + (1-l)*(prob * n2)*((E3[0][i]*S3[0][i])/people) + (1-l)*(prob * n2)*((I3[0][i]*S3[0][i])/people) - alpha*E3[0][i]
    I3[0][i+1] = I3[0][i] + alpha*E3[0][i] - gamma*I3[0][i] - d*I3[0][i]
    R3[0][i+1] = R3[0][i] + gamma*I3[0][i]
    D3[0][i+1] = D3[0][i] + d*I3[0][i]
for i in range(46+12-1,365-1):
    S3[0][i+1] = S3[0][i] - (1-l)*(prob * n2)*((E3[0][i]*S3[0][i])/people) - (1-l)*(prob * n2)*((I3[0][i]*S3[0][i])/people)
    E3[0][i+1] = E3[0][i] + (1-l)*(prob * n2)*((E3[0][i]*S3[0][i])/people) + (1-l)*(prob * n2)*((I3[0][i]*S3[0][i])/people) - alpha*E3[0][i]
    I3[0][i+1] = I3[0][i] + alpha*E3[0][i] - k*gamma*I3[0][i] - d*I3[0][i]
    R3[0][i+1] = R3[0][i] + k*gamma*I3[0][i]
    D3[0][i+1] = D3[0][i] + d*I3[0][i]
max_I3 = np.max(I3[0,:])
with open('各因素叠加的探讨结果.txt','a+') as f:
    f.write('调节社交距离并引入医疗因素和自我防护的结果：\n')
    f.write('最大值：%d\n'%round(max_I3))
    for i in range(len(I3[0,:])):
        if I3[0][i] == max_I3:
            index_max = i
    for i in range(0,365):
        if I3[0][i] <= control_limit and i > index_max and E3[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n'%index)
    r = (round(max_I2) - round(max_I3))/round(max_I2)*100
    f.write('该措施的提升能力(百分比)：%.2f\n\n'%r)
f.close()
#社交距离+医疗举措+自我防护+隔离
S4 = np.zeros((1,365))
E4 = np.zeros((1,365))
I4 = np.zeros((1,365))
R4 = np.zeros((1,365))
D4 = np.zeros((1,365))
S4[0][0] = people - 1
E4[0][0] = 0
I4[0][0] = 1
R4[0][0] = 0
D4[0][0] = 0
for i in range(0,46-1):
    f = I4[0][i]/M 
    if f <= 1:
        w = 2
    else:
        w = 2+n1*(f-1)
    S4[0][i+1] = S4[0][i] - (prob * n1)*((E4[0][i]*S4[0][i])/people) - (prob * w)*((I4[0][i]*S4[0][i])/people)
    E4[0][i+1] = E4[0][i] + (prob * n1)*((E4[0][i]*S4[0][i])/people) + (prob * w)*((I4[0][i]*S4[0][i])/people) - alpha*E4[0][i]
    I4[0][i+1] = I4[0][i] + alpha*E4[0][i] - gamma*I4[0][i] - d*I4[0][i]
    R4[0][i+1] = R4[0][i] + gamma*I4[0][i]
    D4[0][i+1] = D4[0][i] + d*I4[0][i]
for i in range(46-1,46+12-1):
    f = I4[0][i]/M
    if f <= 1:
        w = 2
    else:
        w = 2+n2*(f-1)
    S4[0][i+1] = S4[0][i] - (1-l)*(prob * n2)*((E4[0][i]*S4[0][i])/people) - (1-l)*(prob * w)*((I4[0][i]*S4[0][i])/people)
    E4[0][i+1] = E4[0][i] + (1-l)*(prob * n2)*((E4[0][i]*S4[0][i])/people) + (1-l)*(prob * w)*((I4[0][i]*S4[0][i])/people) - alpha*E4[0][i]
    I4[0][i+1] = I4[0][i] + alpha*E4[0][i] - gamma*I4[0][i] - d*I4[0][i]
    R4[0][i+1] = R4[0][i] + gamma*I4[0][i]
    D4[0][i+1] = D4[0][i] + d*I4[0][i]
for i in range(46+12-1,365-1):
    f = I4[0][i]/(M+M1)
    if f <= 1:
        w = 1
    else:
        w = 1+n2*(f-1)
    S4[0][i+1] = S4[0][i] - (1-l)*(prob * n2)*((E4[0][i]*S4[0][i])/people) - (1-l)*(prob * w)*((I4[0][i]*S4[0][i])/people)
    E4[0][i+1] = E4[0][i] + (1-l)*(prob * n2)*((E4[0][i]*S4[0][i])/people) + (1-l)*(prob * w)*((I4[0][i]*S4[0][i])/people) - alpha*E4[0][i]
    I4[0][i+1] = I4[0][i] + alpha*E4[0][i] - k*gamma*I4[0][i] - d*I4[0][i]
    R4[0][i+1] = R4[0][i] + k*gamma*I4[0][i]
    D4[0][i+1] = D4[0][i] + d*I4[0][i]
max_I4 = np.max(I4[0,:])
with open('各因素叠加的探讨结果.txt','a+') as f:
    f.write('调节社交距离并引入医疗因素和自我防护及隔离的结果：\n')
    f.write('最大值：%d\n'%round(max_I4))
    for i in range(len(I4[0,:])):
        if I4[0][i] == max_I4:
            index_max = i
    for i in range(0,365):
        if I4[0][i] <= control_limit and i > index_max and E4[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n'%index)
    r = (round(max_I3) - round(max_I4))/round(max_I3)*100
    f.write('该措施的提升能力(百分比)：%.2f\n\n'%r)
f.close()
#再次加入小区管控
S5 = np.zeros((1,365))
E5 = np.zeros((1,365))
I5 = np.zeros((1,365))
R5 = np.zeros((1,365))
D5 = np.zeros((1,365))
S5[0][0] = people - 1
E5[0][0] = 0
I5[0][0] = 1
R5[0][0] = 0
D5[0][0] = 0
for i in range(0,46-1):
    f = I5[0][i]/M 
    if f <= 1:
        w = 2
    else:
        w = 2+n1*(f-1)
    S5[0][i+1] = S5[0][i] - (prob * n1)*((E5[0][i]*S5[0][i])/people) - (prob * w)*((I5[0][i]*S5[0][i])/people)
    E5[0][i+1] = E5[0][i] + (prob * n1)*((E5[0][i]*S5[0][i])/people) + (prob * w)*((I5[0][i]*S5[0][i])/people) - alpha*E5[0][i]
    I5[0][i+1] = I5[0][i] + alpha*E5[0][i] - gamma*I5[0][i] - d*I5[0][i]
    R5[0][i+1] = R5[0][i] + gamma*I5[0][i]
    D5[0][i+1] = D5[0][i] + d*I5[0][i]
for i in range(46-1,46+12-1):
    f = I5[0][i]/M
    if f <= 1:
        w = 2
    else:
        w = 2+n2*(f-1)
    S5[0][i+1] = S5[0][i] - (1-l)*(prob * n2)*((E5[0][i]*S5[0][i])/people) - (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people)
    E5[0][i+1] = E5[0][i] + (1-l)*(prob * n2)*((E5[0][i]*S5[0][i])/people) + (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people) - alpha*E5[0][i]
    I5[0][i+1] = I5[0][i] + alpha*E5[0][i] - gamma*I5[0][i] - d*I5[0][i]
    R5[0][i+1] = R5[0][i] + gamma*I5[0][i]
    D5[0][i+1] = D5[0][i] + d*I5[0][i]
for i in range(46+12-1,46+12+7-1):
    f = I5[0][i]/(M+M1)
    if f <= 1:
        w = 1
    else:
        w = 1+n2*(f-1)
    S5[0][i+1] = S5[0][i] - (1-l)*(prob * n2)*((E5[0][i]*S5[0][i])/people) - (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people)
    E5[0][i+1] = E5[0][i] + (1-l)*(prob * n2)*((E5[0][i]*S5[0][i])/people) + (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people) - alpha*E5[0][i]
    I5[0][i+1] = I5[0][i] + alpha*E5[0][i] - k*gamma*I5[0][i] - d*I5[0][i]
    R5[0][i+1] = R5[0][i] + k*gamma*I5[0][i]
    D5[0][i+1] = D5[0][i] + d*I5[0][i]
for i in range(46+12+7-1,365-1):#小区管控
    f = I5[0][i]/(M+M1)
    if f <= 1:
        w = 1
    else:
        w = 1+n3*(f-1)
    S5[0][i+1] = S5[0][i] - (1-l)*(prob * n3)*((E5[0][i]*S5[0][i])/people) - (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people)
    E5[0][i+1] = E5[0][i] + (1-l)*(prob * n3)*((E5[0][i]*S5[0][i])/people) + (1-l)*(prob * w)*((I5[0][i]*S5[0][i])/people) - alpha*E5[0][i]
    I5[0][i+1] = I5[0][i] + alpha*E5[0][i] - k*gamma*I5[0][i] - d*I5[0][i]
    R5[0][i+1] = R5[0][i] + k*gamma*I5[0][i]
    D5[0][i+1] = D5[0][i] + d*I5[0][i]
max_I5 = np.max(I5[0,:])
with open('各因素叠加的探讨结果.txt','a+') as f:
    f.write('进一步引入小区管控的结果：\n')
    f.write('最大值：%d\n'%round(max_I5))
    for i in range(len(I5[0,:])):
        if I5[0][i] == max_I5:
            index_max = i
    for i in range(0,365):
        if I5[0][i] <= control_limit and i > index_max and E5[0][i] <= control_limit:
            index = i+1
            break
    f.write('基本控制住的时间/天：%d\n'%index)
    r = (round(max_I4) - round(max_I5))/round(max_I4)*100
    f.write('该措施的提升能力(百分比)：%.2f\n\n'%r)
f.close()
time = np.array([i for i in range(0,365)]).reshape(1,365)
plt.figure(1,figsize=(8, 6)) 
plt.plot(time[0,:],I0[0,:],color='skyblue')
plt.plot(time[0,:],I1[0,:],color='orange')
plt.plot(time[0,:],I2[0,:],color='hotpink')
plt.plot(time[0,:],I3[0,:],color='red')
plt.plot(time[0,:],I4[0,:],color='green')
plt.plot(time[0,:],I5[0,:],color='blue')
plt.legend(['自然条件','社交距离','社交距离+医疗','社交+医疗+自我防护','社交+医疗+自我+隔离','社交+医疗+自我+隔离+小区管控'])
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('各因素对于患者数量的影响.png')

plt.figure(2,figsize=(8, 6)) 
plt.plot(time[0,:],I4[0,:],color='green',linewidth=0.5)
plt.plot(time[0,:],I5[0,:],color='blue',linewidth=0.5)
plt.legend(['社交+医疗+自我+隔离','社交+医疗+自我+隔离+小区管控'])
plt.xlabel('距离2019.12.8的时间')
plt.ylabel('人数/人')
plt.savefig('是否引入小区管控对于患者数量的影响.png')