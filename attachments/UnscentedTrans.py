
import numpy as np
import matplotlib.pyplot as plt
import random

omega=1/4
r=np.sqrt(2)


#参数初始化--------------------
rho=1
theta=1

rms_rho=0.2
rms_theta=0.8




C=np.array([ [rms_rho*rms_rho, 0], [0 , rms_theta*rms_theta] ])



#无迹变换 (Unscented transform)----------------------------------------
rhotheta=np.zeros(shape=(4,2)) #变换前sigma点

rhotheta[0]=[rho,theta] + r*np.array([rms_rho,0])
rhotheta[1]=[rho,theta] - r*np.array([rms_rho,0])
rhotheta[2]=[rho,theta] + r*np.array([0,rms_theta])
rhotheta[3]=[rho,theta] - r*np.array([0,rms_theta])

xy=np.zeros(shape=(4,2)) #变换后sigma点

for i in range(4):
    xy[i]=[rhotheta[i,0]*np.cos(rhotheta[i,1]), rhotheta[i,0]*np.sin(rhotheta[i,1])]


mean_xy=xy.sum(axis=0) * omega #均值
var_x=0
var_y=0
var_xy=0
for i in range(4): #方差和协方差
    var_x += omega * (xy[i,0]-mean_xy[0])**2
    var_y += omega * (xy[i,1]-mean_xy[1])**2
    var_xy+= omega * (xy[i,0]-mean_xy[0])*(xy[i,1]-mean_xy[1])


#线性变换-------------------------------------------------------------
mean_xl=rho*np.cos(theta)
mean_yl=rho*np.sin(theta)

#变换矩阵
U=np.array([ [np.cos(theta), -rho* np.sin(theta)], [np.sin(theta) , rho* np.cos(theta)] ])
#协方差矩阵
C_l=np.matmul(U,np.matmul(C,U.T))


#蒙特卡洛模拟------------------------------------------------------------
N_mc=100000
rho_mc=np.zeros(N_mc)
theta_mc=np.zeros(N_mc)
for i in range(N_mc):
    rho_mc[i]=random.gauss(rho,rms_rho)
    theta_mc[i]=random.gauss(theta,rms_theta)

x_mc=rho_mc*np.cos(theta_mc)
y_mc=rho_mc*np.sin(theta_mc)

mean_x_mc=x_mc.mean()
mean_y_mc=y_mc.mean()

var_x_mc = ((x_mc-mean_x_mc)**2).mean()
var_y_mc = ((y_mc-mean_y_mc)**2).mean()
var_xy_mc = ((x_mc-mean_x_mc)*(y_mc-mean_y_mc)).mean()



#作图
fig=plt.figure(figsize=(5,5) )
ax=fig.add_subplot(1,1,1)
plt.plot(x_mc,y_mc,'b.',markersize=2  )
plt.xlim(-1.5,1.5)
plt.ylim(-1.5,1.5)

plt.plot(xy[:,0],xy[:,1],'r.' ,markersize=20 )
plt.plot(0,0,'k.' ,markersize=20 )


print('线性变换')
print(f'{mean_xl:.6f}   {mean_yl:.6f}')
print(f'   {C_l[0,0]:.6f}   {C_l[1,1]:.6f}   {C_l[0,1]:.6f}')
print('无迹变换')
print(f'{mean_xy[0]:.6f}   {mean_xy[1]:.6f}')
print(f'   {var_x:.6f}   {var_y:.6f}   {var_xy:.6f}')
print('蒙特卡洛模拟')
print(f'{mean_x_mc:.6f}   {mean_y_mc:.6f}')
print(f'   {var_x_mc:.6f}   {var_y_mc:.6f}   {var_xy_mc:.6f}')



