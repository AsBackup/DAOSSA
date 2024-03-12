
import numpy as np
from astropy.time import Time
from astropy.coordinates import Longitude, Latitude, Angle
from astropy import coordinates
import astropy.units as u


def cart2spher(x,y,z):
    """直角坐标到球坐标（经度、纬度、距离）转换"""
    rho=np.sqrt(x**2+y**2+z**2)
    lon=Longitude(np.arctan2(y/rho,x/rho))
    lat=Latitude(np.arcsin(z/rho))
    return rho,lon,lat

#时间
ot=Time(59654.4625000,format='mjd',scale='utc')

#测站坐标
site_lon=118.78*u.deg
site_lat=32.07*u.deg
site_coord=coordinates.EarthLocation.from_geodetic(site_lon,site_lat,0*u.m)  # 由 测地经纬度 获得 地固系XYZ坐标，同书上(1.64)式


#地心瞬时真赤道真春分点（True Equator and True Equinox，TETE）坐标
r_tete=np.array([-915.731, 5958.096, 3050.295])*u.km



# =================== 基于 ITRS 转换

#恒星时
sg=ot.sidereal_time('apparent','greenwich')

# 由 TETE 转换到 ITRS
r_itrs = np.zeros(3)*u.km # 地心ITRS XYZ坐标（初始化）
cossg=np.cos(sg)
sinsg=np.sin(sg)
r_itrs[0] = r_tete[0] * cossg    + r_tete[1] * sinsg
r_itrs[1] = r_tete[0] * (-sinsg) + r_tete[1] * cossg
r_itrs[2] = r_tete[2]

#站心 ITRS 坐标
sat_topoitrs_xyz = r_itrs - site_coord.get_itrs().cartesian.xyz


#从 站心ITRS 转换到 地平坐标系，基于书(1.53)式变化
cos1 = np.cos(site_lon - 180 * u.deg).value
sin1 = np.sin(site_lon - 180 * u.deg).value
cos2 = np.cos(site_lat - 90 * u.deg).value
sin2 = np.sin(site_lat - 90 * u.deg).value
sat_ah_x = cos2 * cos1 * sat_topoitrs_xyz[0] + cos2 * sin1 * sat_topoitrs_xyz[1] - sin2 * sat_topoitrs_xyz[2]
sat_ah_y = -(-sin1 * sat_topoitrs_xyz[0] + cos1 * sat_topoitrs_xyz[1])
sat_ah_z = sin2 * cos1 * sat_topoitrs_xyz[0] + sin2 * sin1 * sat_topoitrs_xyz[1] + cos2 * sat_topoitrs_xyz[2]

# 由 地平坐标系 XYZ 计算 方位角 azimuth 和 高度角 altitude
sats_rho, sats_az, sats_alt = cart2spher(sat_ah_x, sat_ah_y, sat_ah_z)

print(sats_az.deg, sats_alt.deg)



# =================== 基于 TETE 转换

# 获得测站 GCRS坐标，然后转化到 TETE坐标
site_tete = site_coord.get_gcrs(obstime=ot).transform_to(coordinates.TETE(obstime=ot))

# 获得目标站心赤道坐标系 XYZ 坐标
sats_topotete_cart = r_tete - site_tete.cartesian.xyz

# 计算目标站心赤道坐标系球坐标（赤经、赤纬）
sats_topotete, sats_topotete_ra, sats_topotete_dec = cart2spher(sats_topotete_cart[0], sats_topotete_cart[1], sats_topotete_cart[2])


# 从 站心赤道坐标系 转换到 地平坐标系，基于书（1.54）式变化
cos_sat_dec = np.cos(sats_topotete_dec)
sin_sat_dec = np.sin(sats_topotete_dec)
cos_site_dec = np.cos(site_tete.dec)
sin_site_dec = np.sin(site_tete.dec)
cos_d_ra = np.cos(site_tete.ra - sats_topotete_ra)
sats_A = np.arctan2( -cos_sat_dec * np.sin( site_tete.ra - sats_topotete_ra ),
                    cos_site_dec * sin_sat_dec - sin_site_dec * cos_sat_dec * cos_d_ra )
sats_h = np.arcsin( sin_site_dec * sin_sat_dec + cos_site_dec * cos_sat_dec *  cos_d_ra  )

print(Longitude(sats_A).to_value('deg'), Latitude(sats_h).to_value('deg') )





print('------以下为附加思考，不计入课程评分------')


# ============使用 astropy 的结果

#方法一：直接转换
sat_tete=coordinates.SkyCoord(frame='tete',x=r_tete[0],y=r_tete[1],z=r_tete[2], \
                              representation_type='cartesian',obstime=ot)

sat_altaz_direct = sat_tete.transform_to( coordinates.AltAz( location=site_coord ) )
print(sat_altaz_direct.az.deg,sat_altaz_direct.alt.deg)



#方法二：通过多个站心坐标系转换
sat_itrs=sat_tete.transform_to('itrs')
#站心ITRS坐标
sat_topoitrs_xyz = sat_itrs.cartesian-site_coord.get_itrs().cartesian

#多个站心坐标系转换
sat_topoitrs=coordinates.SkyCoord(sat_topoitrs_xyz,frame='itrs', \
                       representation_type='cartesian',obstime=ot)
sat_cirs=sat_topoitrs.transform_to('cirs')
sat_topocirs=coordinates.SkyCoord(sat_cirs.cartesian,frame='cirs', \
                       representation_type='cartesian',obstime=ot,location=site_coord)

#转换到地平坐标系
sat_topocirs.representation_type='spherical'
sat_altaz=sat_topocirs.transform_to( coordinates.AltAz( obstime=ot, location=site_coord ) )

print(sat_altaz.az.deg,sat_altaz.alt.deg)


#   单纯由于软件包的设计问题，近地目标不能直接在ITRS中加入location后转到altaz
#   详见 https://docs.astropy.org/en/stable/coordinates/common_errors.html
