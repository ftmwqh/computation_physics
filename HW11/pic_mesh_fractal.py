import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename,L):
    f=open(filename,'r')
    lines = f.readline()


    list = lines.strip('\n').split('           ')
    list.remove('')
    x=np.arange(1,L+1,1)
    y=np.arange(1,L+1,1)
    x,y=np.meshgrid(x,y)
    I=np.zeros((L,L))
    for i in range(0,L-1):
        I[i]=list[i*L:(i+1)*L]
        # 把文件读入到一个矩阵
    print(I)
    plt.imshow(I,cmap=plt.cm.gray)
    plt.show()

pic("mesh_2.txt",1001)
pic("mesh_3.txt",1001)
pic("mesh_4.txt",1001)
pic("mesh_5.txt",1001)

