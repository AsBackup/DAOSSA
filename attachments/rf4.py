
import numpy as np
import matplotlib.pyplot as plt

def rk4(x, t, h, f):
    """ h: 步长
        f: 微分方程
    """
    k1 = h * f(x, t)
    k2 = h * f(x + 0.5 * k1, t + 0.5 * h)
    k3 = h * f(x + 0.5 * k2, t + 0.5 * h)
    k4 = h * f(x + k3, t + h)
    return x + (k1 + 2 * k2 + 2 * k3 + k4) / 6


def orbit_motion_equation(x,t):
    """轨道运动方程（无量纲化）
        x[0:3]: 坐标
        x[3:6]: 速度
    """
    r = np.linalg.norm(x[:3])
    y = np.zeros(6) #初始化
    y[:3] = x[3:]
    y[3:] = -x[:3] / r**3
    return y


if __name__=='__main__':

    # 初值
    t0 = 0
    x = np.array([-0.9, -0.4, 0.1, 0.3, -0.9, -0.2])
    h = 0.001

    # 输出节点
    ts = np.linspace(0,10,101)
    ys = np.zeros(shape=(len(ts),6))

    # 数值积分
    t = t0
    for i in range(10100):

        # 搜索输出节点并记录
        t_node = np.where( np.abs(ts - t)<1e-10 )[0]
        if len(t_node)>0:
            ys[t_node[0]] = x.copy()

        x = rk4(x, t, h, orbit_motion_equation)
        t += h


    # 作图
    fig = plt.figure()
    plt.plot(ts, ys, label='runge_kutta')
    plt.show

    fig = plt.figure()
    fig.add_subplot(projection='3d')
    plt.plot(ys[:,0],ys[:,1],ys[:,2], '.-')
    plt.show




