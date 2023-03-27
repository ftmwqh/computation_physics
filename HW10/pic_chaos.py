import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename,a,b,delta):
    f=open(filename,'r')
    lines = f.readlines()
    n=(b-a)/delta+1#总的lambda1取值点个数
    m = len(lines)-int(n)#减去逗号数量
    print('height=', m)

    points=np.zeros((2,m),dtype=float)#记录每个点坐标(lambda,x)
    l=a
    i = 0
    for line in lines:
        list = line.strip('\n') # 除去换行符
        if list==' ,':
            l+=delta
            continue
        points[0][i]=l
        points[1][i]=list
        i += 1
        # 把文件读入到一个矩阵

    fig = plt.figure()
    ax1 = fig.add_subplot(111)


    ax1.scatter(points[0],points[1],s=0.5)
    ax1.set_xlabel(r'$\lambda$')
    ax1.set_ylabel('y')

    plt.show()

pic('result1.txt',-10,10,0.01)
pic('result2.txt',0,2,0.0001)
pic('result3.txt',-2,0,0.0001)