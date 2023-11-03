
"""
基于 astropy.units 的无量纲化示例
更多可参考：https://docs.astropy.org/en/stable/units/index.html
"""

import numpy as np
from astropy import units as u

# JGM-3
u_L = 6378.1363 * u.km
GE = 398600.44150 * u.km**3 / u.s**2
u_T = np.sqrt( u_L**3 / GE )

# 一条轨道根数，所有变量都对单位进行处理
sa = 6800*u.km / u_L
se = 0.1
si = np.pi/4*u.rad
somega = 0.5*u.rad
sw = 5000*u.arcsec
scm = 20*u.deg

# 平运动速度
sn = np.sqrt( 1 / sa**3 ) *u.rad

# numpy 兼容 units 的角度单位，因此进行角度计算时无需转换单位
a = np.pi*u.rad
b = 180*u.deg
print( '三角函数:',np.sin( a ) == np.sin( b ) )

# units 的角度量可以直接加减
slambda = somega + sw + scm

# 同类单位可以直接乘除
scm_t = scm + sn * 10 *u.s / u_T # 10秒后的平近点角
r = 7000000 * u.m / u_L  # 距离
v = np.sqrt ( 2 / r - 1 / sa )  # 活力公式

# 两个 units 相除后不带单位的变量仍然保留 units 属性
x = 5000 * u.km / u_L
y = 4000 * u.km / u_L
alpha = np.arccos( y / x )
print( 'alpha:', alpha ) # 可以看到是自带 u.rad 单位的，因此在后续计算中不需要特别针对反三角函数的结果进行单位转换。

# 输出
print( '输出展示' )
print( v ) # 如果不用输出可以不需要理会它的值和单位
print( v.decompose() ) # 输出时单位化简
print( v.decompose() * u_L / u_T ) # 输出标准单位

# 仅输出值而不带 units 单位时建议使用 .to_value() 而不是 .value ，后者可能因为单位不确定而造成结果错误
print( '不带单位结果\n',slambda.value, slambda.to_value('rad'), slambda.to_value('deg') )



