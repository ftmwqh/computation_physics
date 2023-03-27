import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
plt.rcParams['font.family']=['sans-serif']
plt.rcParams['font.sans-serif']=['Microsoft YaHei']
def pic(filename):
    f=open(filename,'r')
    lines = f.readlines()
    m = len(lines)
    print('height=', m)
    n = len(lines[1].split(' '))
    print('width=', n)
    x = np.zeros(m, dtype=float)
    y=np.zeros(m,dtype=float)
    z=np.zeros(m,dtype=float)
    i = 0
    for line in lines:
        list = line.strip('\n').split(' ')#fortran写入文件，同一行两个数之间是这么多个空格
        try:

            list.remove(' ')
        except ValueError:
            pass

        x[i] = list[0]
        y[i]=list[1]
        z[i]=list[2]
        i += 1
        # 把文件读入到一个矩阵



    fig=plt.figure()
    ax=Axes3D(fig)
    ax.scatter(x,y,z)
    ax.get_proj=lambda :np.dot(Axes3D.get_proj(ax),np.diag([1,1,1.2,1]))#由于作图默认下会有一定坐标轴比例拉伸，这里进行修改
    plt.show()

print(np.sin(np.pi/2))
pic("xy.txt")