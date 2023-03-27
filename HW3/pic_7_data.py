import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']

def pic(filename_1,filename_2,filename_3):
    f=open(filename_1,'r')
    lines = f.readlines()
    m = len(lines)
    print('height=', m)
    n = len(lines[1].split('\t'))
    print('width=', n)
    data_x = np.zeros(m-1, dtype=float)
    data_y=np.zeros(m-1,dtype=float)


    i = 0
    for line in lines:
        if i!=0:
            list = line.strip('\n').split('\t')
            data_x[i-1] = list[0]
            data_y[i-1]=list[1]
        i += 1
        # 把文件读入到一个矩阵
    num_point = sum(data_y)
    data_y = data_y / num_point#归一化

    f_2 = open(filename_2, 'r')
    lines = f_2.readlines()
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
    min_x = min(data)
    max_x = max(data)
    x = np.arange(min_x, max_x, 1)+1  # 找到随机点的最小，最大值，形成分布区间
    mine = np.zeros(int(max_x-min_x))  # 统计随机点分布

    for i in range(0, len(data)):
        mine[int(data[i]-min_x-1)]+=1
    mine=mine/len(data)
    num_point=sum(data_y)
    data_y=data_y/num_point

    f_3 = open(filename_3, 'r')#读取直接抽样法数据
    lines = f_3.readlines()
    m = len(lines)
    print('height=', m)
    n = len(lines[1].split(' '))
    print('width=', n)
    data_2 = np.zeros(m, dtype=float)

    i = 0
    for line in lines:
        list = line.strip('\n')
        data_2[i] = list
        i += 1
        # 把文件读入到一个矩阵
    print(data_2)
    min_x = min(data_2)
    max_x = max(data_2)
    x_2 = np.arange(min_x, max_x, 1) +1  # 找到随机点的最小，最大值，形成分布区间
    print(x_2)
    straight = np.zeros(int(max_x - min_x))  # 统计随机点分布

    for i in range(0, len(data_2)):
        straight[int(data_2[i] - min_x - 1)] += 1
    straight = straight / len(data_2)


    plt.plot(x, mine, label='舍选法',linewidth=1)
    plt.plot(data_x,data_y, label='标准结果',linewidth=1)
    plt.plot(x_2,straight,label='直接抽样法',linewidth=1)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend(loc='upper left')
    plt.show()

pic("data.TXT","x_1.txt","x_2.txt")