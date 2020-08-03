#### **一、模型假设：**

1. 武汉市在疫情防控期间没有外来流动人口，只考虑本市的常住人口
2. 新冠肺炎病毒的自然传播能力是一定的
3. 支援武汉的医疗力量是国家规定的固定值，不能改变
4. 支援武汉的医疗力量只考虑病床的增加，不考虑医务人员参与对于疫情的影响
5. 所有增加的病床均在2020.2.4投入使用，忽略真实情况的火神山、雷神山和方舱医院不同的开启时间
6. 新冠肺炎只有在密切接触时才可以传播
7. 武汉市的人口均匀分布
8. 最严格的管控条件下，密切接触人数为4
9. 新冠肺炎患者的康复期为14天
10. 新冠肺炎的潜伏期为7天
11. SIRD和SEIRD模型中I的最大值对应于现实生活中的累计确诊数量
12. 武汉市解封时，疫情基本得到控制

#### **二、变量汇总：**

<img src="https://github.com/ziseSunny/virus/blob/master/picture/变量说明1.png" alt="变量说明1" style="zoom:50%;" >
<img src="https://github.com/ziseSunny/virus/blob/master/picture/变量说明2.png" alt="变量说明2" style="zoom:50%;" >

#### **三、变量取值：**

<img src="https://github.com/ziseSunny/virus/blob/master/picture/变量取值.png" alt="变量取值" style="zoom:50%;" >
#### **四、模型构建：**

##### **（一）SIRD模型（在SIR的基础上单独考虑D）：**

1. 示意图：

   **SIRD模型的示意图：**

   <img src="https://github.com/ziseSunny/virus/blob/master/封城时间对于患者数量的影响.png" alt="SIR模型" style="zoom:30%;" >

   **密切接触的模型：**
   整个正方形代表所在区域，圆形代表人

   <img src="https://github.com/ziseSunny/virus/blob/master/密切接触.png" alt="密切接触" style="zoom:10%;" />

2. **密切接触的公式：**
   n的公式通过投点实验归纳得出
   <img src="https://github.com/ziseSunny/virus/blob/master/picture/密切接触公式.png" alt="密切接触公式" style="zoom:30%;" >

   每个圆形都具有流动性，故n包含了中心点
   
3. **1阶段模型：**

   <img src="https://github.com/ziseSunny/virus/blob/master/picture/1阶段模型.png" alt="1阶段模型" style="zoom:50%;" >

4. **2阶段模型：**

   <img src="https://github.com/ziseSunny/virus/blob/master/picture/2阶段模型.png" alt="2阶段模型" style="zoom:50%;" >
5. **3阶段模型：**

   <img src="https://github.com/ziseSunny/virus/blob/master/picture/变量说明1.png" alt="3阶段模型" style="zoom:50%;" >

**（二）SEIRD模型（在SEIR的基础上单独考虑D）：**

1. 参数值的处理：

   根据$SIRD$模型，得到第1阶段和第2阶段对应的p，第$2$阶段相对较崩溃，第1阶段还处于疫情防控的初期，故此时采用的p为两者的平均值，即p = (0.01671555+0.05222604)/2 = 0.034470795$

   在未引入流动系数时，第2阶段和第3阶段的n=4

2. 示意图：$\beta$代表的是易感人群密切接触潜伏期人群或患者而成为潜伏期患者的概率

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/SEIRD.png" alt="SEIRD" style="zoom:50%;" />

3. **基本的$3$阶段模型：**

   此模型没有引入自我防护、隔离、流动系数和小区管控

   **微分形式：**

   第$1$阶段：
   $$
   n = \lceil \frac {d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} \frac{dS}{dt} = -\beta \frac{SE}{N} -\beta \frac{SI}{N} \\ \frac{dE}{dt} = \beta\frac{SE}{N} + \beta\frac{SI}{N} - \alpha E \\ \frac{dI}{dt} = \alpha E - \gamma I - dI \\ \frac{dR}{dt} = \gamma I \\ \frac{dD}{dt} = dI  \end{cases}
   $$

   第$2$阶段：
   $$
   n = 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} \frac{dS}{dt} = -\beta \frac{SE}{N} -\beta \frac{SI}{N} \\ \frac{dE}{dt} = \beta\frac{SE}{N} + \beta\frac{SI}{N} - \alpha E \\ \frac{dI}{dt} = \alpha E - \gamma I - dI \\ \frac{dR}{dt} = \gamma I \\ \frac{dD}{dt} = dI  \end{cases}
   $$

   第$3$阶段：
   $$
   n = 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} \frac{dS}{dt} = -\beta \frac{SE}{N} -\beta \frac{SI}{N} \\ \frac{dE}{dt} = \beta\frac{SE}{N} + \beta\frac{SI}{N} - \alpha E \\ \frac{dI}{dt} = \alpha E - k\gamma I - dI \\ \frac{dR}{dt} = k\gamma I \\ \frac{dD}{dt} = dI  \end{cases}
   $$

   由于$p$、$n$、$\alpha$、$\gamma$和$k$均已知，直接采用差分形式进行模拟

   **差分形式：**

   第$1$阶段：
   $$
   n = \lceil \frac {d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + \beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$2$阶段：
   $$
   n = 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + \beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$3$阶段：
   $$
   n = 4
   $$

   $$
   \beta = pn
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + \beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - k\gamma I_i - dI_i \\ R_{i+1} = R_i + k\gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   判断疫情基本得到控制的条件为：$I \leq c$ && $E \leq c$ && 此时的时间大于峰值对应的时间

4. **加入自我防护的$3$阶段模型：**

   自我防护措施主要是引入变量$l$，综合考虑口罩等措施的实际防护能力和物资紧缺的状况，令$l = 0.7$

   假设封城开始后，人们就自觉进行自我防护

   **差分形式：**

   第$1$阶段：
   $$
   n = \lceil \frac {d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = np
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + \beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$2$阶段：
   $$
   n = 4
   $$

   $$
   \beta = np
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)\beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$3$阶段：
   $$
   n = 4
   $$

   $$
   \beta = np
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)\beta \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)\beta\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - k\gamma I_i - dI_i \\ R_{i+1} = R_i + k\gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

5. **加入隔离的$3$阶段模型：**

   隔离主要影响易感者密切接触的患者数，此时用$w$来衡量易感者在不同阶段接触的患者数量

   与$w$有关的变量$f$的表达式为：
   $$
   f = \frac {I}{M+M'}
   $$
   **差分形式：**

   第$1$阶段：
   $$
   n = \lceil \frac {d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = np
   $$

   $$
   w = \begin{cases} w_0 & f \leq 1 \\ 
                     w_0+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$2$阶段：
   $$
   n = 4
   $$

   $$
   \beta = np
   $$

   $$
   w = \begin{cases} w_0 & f \leq 1 \\ 
                     w_0+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$3$阶段：
   $$
   n = 4
   $$

   $$
   \beta = np
   $$

   $$
   w = \begin{cases} w_1 & f \leq 1 \\ 
                     w_1+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - k\gamma I_i - dI_i \\ R_{i+1} = R_i + k\gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

6. **流动系数的处理：**

   流动系数主要影响密切接触的距离，当所在区域的流动性增强时，密切接触的距离会增大

   考虑在封城开始后，所在区域的流动系数发生变化

   流动系数$t$的表达式为：
   $$
   t = \begin{cases} 1 & 封城前 \\ 
                     0.3 & 封城后 \end{cases}
   $$

7. **加入小区管控的$4$阶段模型：**

   根据相关资料，武汉市从$2020.2.11$开始全面实施小区管控，对应的$4$个阶段为：

   $2019.12.8～2020.1.22$：共$46$天

   $2020.1.23～2020.2.3$：共$12$天

   $2020.2.4～2020.2.10$：共$7$天

   $2020.2.11～$

   **差分形式：**

   第$1$阶段：
   $$
   t = 1
   $$

   $$
   n = \lceil \frac {t * d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = np
   $$

   $$
   w = \begin{cases} w_0 & f \leq 1 \\ 
                     w_0+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -\beta \frac{S_iE_i}{N} -wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + \beta\frac{S_iE_i}{N} + wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$2$阶段：
   $$
   t = 0.3
   $$

   $$
   n = \lceil \frac {t * d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   \beta = np
   $$

   $$
   w = \begin{cases} w_0 & f \leq 1 \\ 
                     w_0+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - \gamma I_i - dI_i \\ R_{i+1} = R_i + \gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$3$阶段：
   $$
   t = 0.3
   $$

   $$
   n = \lceil \frac {t * d_t}{2\rho} \rceil^2 + 4
   $$

   $$
   w = \begin{cases} w_1 & f \leq 1 \\ 
                     w_1+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - k\gamma I_i - dI_i \\ R_{i+1} = R_i + k\gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   第$4$阶段：
   $$
   n = 4
   $$

   $$
   w = \begin{cases} w_1 & f \leq 1 \\ 
                     w_1+n(f-1) & f > 1 \end{cases}
   $$

   $$
   \beta = np
   $$

   $$
   \begin{cases} S_{i+1} = S_i -(1-l)\beta \frac{S_iE_i}{N} -(1-l)wp \frac{S_iI_i}{N} \\ E_{i+1} = E_i + (1-l)\beta\frac{S_iE_i}{N} + (1-l)wp\frac{S_iI_i}{N} - \alpha E_i \\ I_{i+1} = I_i + \alpha E_i - k\gamma I_i - dI_i \\ R_{i+1} = R_i + k\gamma I_i \\ D_{i+1} = D_i + dI_i  \end{cases}
   $$

   模拟结果为：

   **患者人数的变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/引入流动系数后的SEIR模型对应的患者人数趋势图.png" alt="引入流动系数后的SEIR模型对应的患者人数趋势图" style="zoom:60%;" />

   累计确诊数量：$65861$，疫情得到基本控制的时间为：第$140$天

   根据公开数据，武汉市的累计确诊数量为：$50340$，共封城$76$天，则疫情基本得到控制的时间为：第$46+76=122$天，考虑到新冠肺炎的潜伏期为$7$天，且初期数据不够透明，本模型模拟的结果较为合理

#### **五、敏感性分析：**

**在引入流动系数和小区管控模型的基础上进行讨论：**

1. **封城时间的敏感性分析：**

   **封城时间与累计确诊数量、基本得到控制的时间点：**

   |      封城时间      | 第$37$天  | 第$42$天  | 第$47$天  | 第$52$天  | 第$57$天  |
   | :----------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |    累计确诊数量    |  $6290$   |  $17398$  |  $65861$  | $118843$  | $2464403$ |
   | 基本得到控制的时间 | 第$103$天 | 第$116$天 | 第$140$天 | 第$171$天 | 第$159$天 |

   **不同封城时间对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/封城时间对于患者数量的影响.png" alt="封城时间对于患者数量的影响" style="zoom:60%;" />

   $d_1=56$对应的$I$过大，选取前$4$条曲线单独进行分析：

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/封城时间对于患者数量的影响_clear.png" alt="封城时间对于患者数量的影响_clear" style="zoom:60%;" />

   **封城时间对于患者人数的影响较为显著**

2. **自我防护开始时间的敏感性分析：**

   **自我防护开始时间与累计确诊数量、基本得到控制的时间点：**

   |  自我防护开始时间  | 第$37$天 | 第$42$天  | 第$47$天  | 第$52$天  | 第$57$天  |
   | :----------------: | :------: | :-------: | :-------: | :-------: | :-------: |
   |    累计确诊数量    |  $3977$  |  $13663$  |  $65861$  | $2228386$ | $2802732$ |
   | 基本得到控制的时间 | 第$97$天 | 第$112$天 | 第$140$天 | 第$171$天 | 第$148$天 |

   **不同自我防护开始时间对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/自我防护时间对于患者数量的影响.png" alt="自我防护时间对于患者数量的影响" style="zoom:60%;" />

   $d_2=51$和$d_2=56$对应的$I$太大，选取前$3$条曲线单独进行分析：

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/自我防护时间对于患者数量的影响_clear.png" alt="自我防护时间对于患者数量的影响_clear" style="zoom:60%;" />

   **自我防护开始时间对于患者人数的影响较为显著**

3. **医疗援助到达时间的敏感性分析：**

   **医疗援助到达时间与累计确诊数量、基本得到控制的时间点：**

   |  医疗援助到达时间  | 第$49$天  | 第$54$天  | 第$59$天  | 第$64$天  | 第$69$天  |
   | :----------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |    累计确诊数量    |  $30258$  |  $47353$  |  $65861$  | $102631$  | $1888262$ |
   | 基本得到控制的时间 | 第$122$天 | 第$129$天 | 第$140$天 | 第$160$天 | 第$198$天 |

   **不同医疗援助到达时间对应的患者人数趋势图：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/医疗援助时间对于患者数量的影响.png" alt="医疗援助时间对于患者数量的影响" style="zoom:60%;" />

   $d_3=68$对应的$I$太大，选取前$4$条曲线进行分析：

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/医疗援助时间对于患者数量的影响_clear.png" alt="医疗援助时间对于患者数量的影响_clear" style="zoom:60%;" />

   **医疗援助到达时间对患者人数的影响不如封城时间和自我防护开始时间显著，但$d_3=63～68$中存在一个关键的时间点，一旦医疗援助到达时间晚于此时间点，防控形势将发生很大变化**

4. **小区管控开始时间的敏感性分析：**

   **小区管控开始时间与累计确诊数量、基本得到控制的时间点：**

   |  小区管控开始时间  | 第$56$天  | 第$61$天  | 第$66$天  | 第$71$天  | 第$76$天  |
   | :----------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |    累计确诊数量    |  $62212$  |  $65861$  |  $65861$  |  $65861$  |  $65861$  |
   | 基本得到控制的时间 | 第$133$天 | 第$139$天 | 第$140$天 | 第$141$天 | 第$141$天 |

   **不同小区管控开始时间对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/小区管控时间对于患者数量的影响.png" alt="小区管控时间对于患者数量的影响" style="zoom:60%;" />

   **小区管控是最为严苛的防控措施，虽然不能很好地控制累计确诊数量，但可以有效控制防控进度**

5. **自我防护效果的敏感性分析：**

   **自我防护效果与累计确诊数量、基本得到控制的时间点：**

   |    自我防护效果    |   $0.6$   |  $0.65$   |   $0.7$   |  $0.75$   |   $0.8$   |
   | :----------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |    累计确诊数量    |  $94726$  |  $75328$  |  $65861$  |  $57681$  |  $50606$  |
   | 基本得到控制的时间 | 第$200$天 | 第$156$天 | 第$140$天 | 第$129$天 | 第$121$天 |

   **不同自我防护效果对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/自我防护的效果对于患者数量的影响.png" alt="自我防护的效果对于患者数量的影响" style="zoom:60%;" />

   **自我防护效果可以有效控制疫情的扩散，也可以控制疫情防控的进度**

6. **医疗援助到达前的隔离强度的敏感性分析：**

   **医疗援助到达前的隔离强度与累计确诊数量、基本得到控制的时间点：**

   | 医疗援助到达前的隔离效果 |   $1.5$   |  $1.75$   |    $2$    |  $2.25$   |   $2.5$   |
   | :----------------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |       累计确诊数量       |  $32151$  |  $45725$  |  $65861$  |  $96271$  | $2052016$ |
   |    基本得到控制的时间    | 第$125$天 | 第$132$天 | 第$140$天 | 第$154$天 | 第$181$天 |

   **不同医疗援助到达前的隔离强度对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/医疗援助前的隔离效果对于患者数量的影响.png" alt="医疗援助前的隔离效果对于患者数量的影响" style="zoom:60%;" />

   $w_0=2.5$对应的$I$太大，选取前$4$条曲线进行分析：

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/医疗援助前的隔离效果对于患者数量的影响_clear.png" alt="医疗援助前的隔离效果对于患者数量的影响_clear" style="zoom:60%;" />

   **医疗援助到达前的隔离强度对于患者人数的影响较为显著，且$w_0=2.25～2.5$中存在一个关键点，超过此点，疫情防控形势将发生很大变化**

7. **医疗援助到达后的隔离强度的敏感性分析：**

   **医疗援助到达后的隔离强度与累计确诊数量、基本得到控制的时间点：**

   | 医疗援助到达后的隔离效果 |   $0.5$   |  $0.75$   |    $1$    |  $1.25$   |   $1.5$   |
   | :----------------------: | :-------: | :-------: | :-------: | :-------: | :-------: |
   |       累计确诊数量       |  $65861$  |  $65861$  |  $65861$  |  $65861$  |  $65861$  |
   |    基本得到控制的时间    | 第$136$天 | 第$138$天 | 第$140$天 | 第$142$天 | 第$144$天 |

   **不同医疗援助到达后的隔离强度对应的患者人数变化趋势：**

   <img src="https://github.com/ziseSunny/virus-simulation/blob/master/医疗援助后的隔离效果对于患者数量的影响.png" alt="医疗援助后的隔离效果对于患者数量的影响" style="zoom:60%;" />

   **医疗援助到达后的隔离强度对于患者人数的影响不显著，但可以控制疫情防控的进度**

#### **六、各因素重要性分析：**

在最基本的$SEIRD$模型的基础上逐步增加因素：

1. 加入社交距离的控制（封城措施）
2. 加入医疗因素
3. 加入医疗援助
4. 加入隔离
5. 加入小区管控

**逐步加入各因素对应的累计确诊数量、基本得到控制的时间：**

|           操作步骤           |  无操作   |  步骤$1$  |  步骤$2$  |  步骤$3$  |  步骤$4$  |  步骤$5$  |
| :--------------------------: | :-------: | :-------: | :-------: | :-------: | :-------: | :-------: |
|         累计确诊数量         | $4730700$ | $3635181$ | $3635181$ | $2703591$ |  $65861$  |  $65861$  |
|      基本得到控制的时间      | 第$186$天 | 第$208$天 | 第$150$天 | 第$152$天 | 第$149$天 | 第$140$天 |
| 较于上一步确诊数量控制的比例 |     -     | $23.16$%  |  $0.00$%  | $25.63$%  | $97.56$%  |  $0.00$%  |

**逐步加入各因素对应的患者人数变化趋势：**

<img src="https://github.com/ziseSunny/virus-simulation/blob/master/各因素对于患者数量的影响.png" alt="各因素对于患者数量的影响" style="zoom:60%;" />

最后的$2$条线几近重合，单独进行分析：

<img src="https://github.com/ziseSunny/virus-simulation/blob/master/是否引入小区管控对于患者数量的影响.png" alt="是否引入小区管控对于患者数量的影响" style="zoom:60%;" />

**有关的思考：**

在疫情防控中，最关键的因素为隔离，控制社交距离和医疗援助也有效控制疫情的进一步扩散

医疗因素和小区管控可以加速疫情防控的进程，提升疫情防控的速度

#### **七、群体免疫的思考：**

在$SEIRD$的最终模型中考察群体免疫的存在性：

**现有易感人群数量的变化趋势：**

<img src="https://github.com/ziseSunny/virus-simulation/blob/master/群体免疫的思考.png" alt="群体免疫的思考" style="zoom:60%;" />

在疫情得到控制后，易感人群的数量有所下降，有一部分人群对新冠肺炎已有自身免疫，说明群体免疫是存在的

#### **八、模型创新点：**

1. 在最基本的传染病模型——$SIR$的基础上，结合本次武汉市新冠肺炎疫情的相关数据，创造性地引入医疗因素、密切接触模型和死亡人数，构建$SRID$模型
2. 在$SRID$模型的基础上，考虑到新冠肺炎存在潜伏期，结合$SEIR$模型，引入流动系数和隔离系数、自我防护效果，构建$SERID$模型
3. 在较为简单的基础模型上创造性加入各因素，模拟效果接近真实情形，具有说服力和可信性

#### **九、疫情防控的相关思考和建议：**

1. 有效控制社交距离，合理设置封城时间是疫情初期防控的关键
2. 及时安排医疗援助是在疫情迅猛发展期有效防控的关键
3. 增加隔离效果，有效控制易感人群密切接触的患者数是整个疫情防控过程的关键
4. 小区管控可以加速疫情防控进程，医疗援助对应的医疗因素变化可以加速疫情防控进程
5. 新冠肺炎与人类会共存很长一部分时间，国家一定要加强公共卫生机制和设施的建设，有效实现疫情防控常态化，不断提升应急医疗的建设，不断提升疫情防控能力

#### **十、参考资料和数据来源：**

1. 武汉市疫情数据：https://github.com/839Studio/Novel-Coronavirus-Updates

2. $SIR$模型：

   https://blog.csdn.net/qq_42677001/article/details/104254761?utm_medium=distribute.pc_relevant_t0.none-task-blog-BlogCommendFromMachineLearnPai2-1.channel_param&depth_1-utm_source=distribute.pc_relevant_t0.none-task-blog-BlogCommendFromMachineLearnPai2-1.channel_param

3. 新冠肺炎死亡率和$SIRD$模型：

   蔡洁,贾浩源,王珂.基于SEIR模型对武汉市新型冠状病毒肺炎疫情发展趋势预测[J].山东医药,2020,60(06):1-4.

4. 武汉市病床数量：

   https://new.qq.com/omn/20200128/20200128A0HE4100.html

   $M=600+5311+6000=11911$

   火神山、雷神山医院的病床数：百度百科

   方舱医院的病床数：https://baijiahao.baidu.com/s?id=1673679869890764276&wfr=spider&for=pc

   $M’=1000+1600+13000$

5. 武汉市封城总时间：https://www.sohu.com/a/386310760_779534

6. 武汉市小区管控时间：http://wh.bendibao.com/news/2020211/107372.shtm

7. 武汉市累计确诊数量：https://voice.baidu.com/act/newpneumonia/newpneumonia/?from=osari_aladin_banner&city=湖北-武汉
