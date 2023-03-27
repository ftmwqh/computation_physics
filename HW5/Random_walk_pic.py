import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename1,filename2,filename3,A,w):#绘制随机行走各平均值与期望对比
    f1=open(filename1,'r')
    f2=open(filename2,'r')
    f3 = open(filename3, 'r')
    lines1 = f1.readlines()
    lines2 = f2.readlines()
    lines3 = f3.readlines()
    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    x = np.zeros(m, dtype=float)
    x_standard=np.zeros(m, dtype=float)
    y_standard = np.zeros(m, dtype=float)
    y=np.zeros(m,dtype=float)
    y=np.zeros(m,dtype=float)
    x2 = np.zeros(m, dtype=float)
    x2_standard = np.zeros(m, dtype=float)
    var = np.zeros(m, dtype=float)
    var_standard=np.zeros(m, dtype=float)
    t=np.arange(1,101,1)

    i = 0#读取x模拟值
    for line in lines1:
        list = line.strip('\n')#除去换行符
        x[i] = list
        i += 1
        # 把文件读入到一个矩阵

    i=0#读取y模拟值
    for line in lines2:
        list = line.strip('\n')#除去换行符
        y[i] = list
        i += 1
        # 把文件读入到一个矩阵

    i = 0#读取x^2模拟值
    for line in lines3:
        list = line.strip('\n')  # 除去换行符
        x2[i] = list
        i += 1
        # 把文件读入到一个矩阵

    # 随机游走方差结果
    var=x2-np.power(x,2)-np.power(y,2)
    # x的期望
    x_standard=A/2*np.sin(w/2*t)*np.sin(w/2*(t+1))/np.sin(w/2)
    # x^2的期望
    x2_standard=t*(1-A**2/8)+A**2/8*(np.sin(w*t)*np.cos(w*(t+1)))/np.sin(w)+A**2/4*(np.sin(w*t/2)*np.sin(w*(t+1)/2)/np.sin(w/2))**2
    #方差期望
    var_standard=t*(1-A**2/8)+A**2/8*(np.sin(w*t)*np.cos(w*(t+1))/np.sin(w))

    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3 = fig.add_subplot(223)
    ax4=fig.add_subplot(224)

    ax1.scatter(t, x, label='模拟结果', color='r')
    ax1.plot(t, x_standard, label='<x(t)>', linewidth=1)
    ax1.set_xlabel('t')
    ax1.set_ylabel('x')
    ax1.set_title("(a)<x(t)>模拟",y=-0.2)
    ax1.legend(loc='upper left')

    ax2.scatter(t, y, label='模拟结果', color='r',s=5)
    ax2.plot(t, y_standard, label='<y(t)>', linewidth=1)
    ax2.set_xlabel('t')
    ax2.set_ylabel('y')
    ax2.set_title("(b)<y(t)>模拟",y=-0.2)
    ax2.set_ylim((-1,1))
    ax2.legend(loc='upper left')

    ax3.scatter(t,x2, label='模拟结果', color='r',s=5)
    ax3.plot(t, x2_standard, label=r'$<r^2(t)>$', linewidth=1)
    ax3.set_xlabel('t')
    ax3.set_ylabel(r'$x^2$')
    ax3.set_title(r'(c)$<r^2(t)>$模拟',y=-0.2)
    ax3.legend(loc='upper left')

    ax4.scatter(t, var, label='模拟结果', color='r', s=5)
    ax4.plot(t, var_standard, label=r'$Var_r(t)$', linewidth=1)
    ax4.set_xlabel('t')
    ax4.set_ylabel(r'$Var(x)$')
    ax4.set_title(r'(d)$Var_r(t)$模拟',y=-0.2)
    ax4.legend(loc='upper left')
    plt.show()

def pic_C(filename1,A,w,vx0_ave,vy0_ave,vx02_ave,vy02_ave,tau):#绘制速度自相关函数曲线
    f1=open(filename1,'r')
    lines1 = f1.readlines()
    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    C = np.zeros(m, dtype=float)
    C_standard=np.zeros(m, dtype=float)
    C_approx=np.zeros(m,dtype=float)
    t=np.arange(0,101,1)
    i = 0#读取x模拟值
    for line in lines1:
        list = line.strip('\n')#除去换行符
        C[i] = list
        i += 1
        # 把文件读入到一个矩阵

    C_standard=(vx02_ave+vy02_ave)*np.exp(-t/tau)/2+A/(2*tau)*vx0_ave/2*(tau*np.sin(w*t)-w*tau**2*(np.cos(w*t)-np.exp(-t/tau)))/(1+w**2*tau**2)
    C_approx=A/(2*tau)*vx0_ave/2*tau*np.sin(w*t)
    fig=plt.figure()
    ax1=fig.add_subplot(111)
    print(C_standard[0])
    ax1.scatter(t[0:100], C, label='模拟结果', color='r')
    ax1.plot(t, C_standard, label='C(t)', linewidth=1)#自相关函数标准函数画图
    #ax1.plot(t, C_approx, label='C(t)近似', linewidth=1)#舍去指数衰减项和二阶量的近似曲线
    ax1.set_xlabel('t')
    ax1.set_ylabel('C')
    ax1.legend(loc='upper right')

    plt.show()

def walk(filename1,filename2):#绘制随机行走坐标变化图
    f1=open(filename1,'r')
    f2=open(filename2,'r')

    lines1 = f1.readlines()
    lines2 = f2.readlines()

    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split(' '))
    print('width=', n)
    x = np.zeros(m, dtype=float)
    y=np.zeros(m,dtype=float)
    t=np.arange(1,101,1)

    i = 0#读取x模拟值
    for line in lines1:
        list = line.strip('\n')#除去换行符
        x[i] = list
        i += 1
        # 把文件读入到一个矩阵

    i=0#读取y模拟值
    for line in lines2:
        list = line.strip('\n')#除去换行符
        y[i] = list
        i += 1
        # 把文件读入到一个矩阵

    fig=plt.figure()
    ax1=fig.add_subplot(221)
    ax2=fig.add_subplot(222)
    ax3 = fig.add_subplot(223)

    ax1.plot(t, x, label='模拟结果', marker='o')
    ax1.set_xlabel('t')
    ax1.set_ylabel('x')
    ax1.set_title("(a)x(t)模拟",y=-0.2)
    ax1.legend(loc='upper left')

    ax2.plot(t, y, label='模拟结果', marker='o')
    ax2.set_xlabel('t')
    ax2.set_ylabel('y')
    ax2.set_title("(b)y(t)模拟",y=-0.2)
    ax2.legend(loc='upper left')

    ax3.plot(x,y, label='模拟结果', marker='o')
    ax3.set_xlabel('x')
    ax3.set_ylabel('y')
    ax3.set_title('(c)粒子轨迹',y=-0.2)
    ax3.legend(loc='upper left')

    plt.show()

def statistic(filename_1,filename_2):#统计粒子分布
    f1=open(filename_1,'r')
    lines1 = f1.readlines()
    f2 = open(filename_2, 'r')
    lines2 = f2.readlines()
    m = len(lines1)
    print('height=', m)
    n = len(lines1[1].split('\t'))
    print('width=', n)
    data_x = np.zeros(m, dtype=float)
    data_y=np.zeros(m,dtype=float)


    i = 0
    for line in lines1:
        list = line.strip('\n').split('\t')
        data_x[i] = list[0]
        i += 1

    i = 0
    for line in lines2:
        list = line.strip('\n').split('\t')
        data_y[i] = list[0]
        i += 1
        # 把文件读入到一个矩阵


    x = np.arange(-50, 50, 1)+1
    px = np.zeros(100)  # 统计x坐标分布
    px_standard=np.zeros(100)
    py = np.zeros(100)  #统计y坐标分布
    py_standard = np.zeros(100)

    for data in data_x:
        px[int(data+50-1)]+=1
    for data in data_y:
        py[int(data+50-1)]+=1
    px=px/m
    sigmax2=100*(1/2-1/8)+1/8*np.sin(0.5*100)*np.cos(101*0.5)/np.sin(0.5)
    x_ave=1/2*np.sin(0.5*100/2)*np.sin(101*0.5/2)/np.sin(0.5/2)
    px_standard=1/np.sqrt(2*np.pi*sigmax2)*np.exp(-np.power(x-x_ave,2)/(2*sigmax2))

    py=py/m
    sigmay2=100/2
    y_ave=0
    py_standard=1/np.sqrt(2*np.pi*sigmay2)*np.exp(-np.power(x-y_ave,2)/(2*sigmay2))

    fig = plt.figure()
    ax1 = fig.add_subplot(211)
    ax2 = fig.add_subplot(212)

    ax1.bar(x, px,1, label='模拟结果')
    ax1.plot(x,px_standard, label='理论分布',color='r',linewidth=1)
    ax1.set_xlabel('x')
    ax1.set_ylabel(r'$p_x(x)$')
    ax1.set_title("(a)x方向分布", y=-0.2)
    ax1.legend(loc='upper left')

    ax2.bar(x, py,1, label='模拟结果')
    ax2.plot(x,py_standard, label='理论分布',color='r',linewidth=1)
    ax2.set_xlabel('y')
    ax2.set_ylabel(r'$p_y(y)$')
    ax2.set_title("(b)y方向分布", y=-0.2)
    ax2.legend(loc='upper left')

    plt.show()

pic_C("C_A1w05_1.txt",1,0.5,0,0,1,1,10**(-7))

pic("A1w1x.txt","A1w1y.txt","A1w1x2.txt",1,1)
pic("A1w01x.txt","A1w01y.txt","A1w01x2.txt",1,0.1)
pic("A1w05x.txt","A1w05y.txt","A1w05x2.txt",1,0.5)

pic("A05w05x.txt","A05w05y.txt","A05w05x2.txt",0.5,0.5)
walk('walkx.txt','walky.txt')
statistic('distribution_x.txt','distribution_y.txt')
pic_C("C_A1w05.txt",1,0.5,1,1,1,1,10**(-7))