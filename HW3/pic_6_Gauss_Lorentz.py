import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename,divide_n):#绘图
    f=open(filename,'r')
    lines = f.readlines()
    m = len(lines)
    print('height=', m)
    n = len(lines[1].split(' '))
    print('width=', n)
    data = np.zeros(m, dtype=float)


    i = 0
    for line in lines:
        list = line.strip('\n')
        data[i] = list
        i += 1
        # 把文件读入到一个矩阵
    print(data)
    min_x=min(data)
    max_x=max(data)
    x=np.linspace(min_x,max_x,divide_n)#找到随机点的最小，最大值，形成分布区间
    print(x)
    mine=np.zeros(divide_n-1)#统计随机点分布
    for j in range(0,divide_n-1):
        for i in range(0,len(data)):
            if x[j]<data[i]<x[j+1]:
                mine[j]+=1
    a=3.0
    standard=np.zeros(divide_n-1)
    standard=np.sqrt(a/np.pi)*np.exp(-a*x**2)#标准p函数
    mine=mine/(len(data)*(x[1]-x[0]))#归一化
    x_draw=x[0:divide_n-1]+(x[-1]-x[0])/(2*divide_n-2)#由于记录的是区间统计数量，实际的分布函数格数比区间端点数少1
    plt.plot(x_draw, mine, label='计算结果')
    plt.plot(x, standard, label='标准结果')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc='upper left')
    plt.show()

a=[1,-4,3]
print(min(a))
pic('x.txt',100)