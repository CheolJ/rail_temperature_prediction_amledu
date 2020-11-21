# Library
import math
import datetime
import numpy as np


# constant variable
orbit_days = 365.256363004  # Earth orbit in days
au = 149598261  # The semi-major axis of the oribital ellipse
e = 0.01671123  # Earth orbit elliptical eccentricity
solar_c = 1367  # Solar constant

def el_az_changer(time_raw, lat = 36, lte = 127):

    hour = time_raw.hour + time_raw.minute / 60
    do = datetime.datetime(time_raw.year, 1, 1)
    d1 = datetime.datetime(time_raw.year, time_raw.month, time_raw.day, time_raw.hour)
    day_delta = d1 - do
    day = day_delta.days

    dy = time_raw.year - 1949
    leap = round(dy / 4)
    Jd = 2432916.5 + dy * 365 + leap + day + hour / 24  # 율리우스 적일
    n = Jd - 2451545.0

    L = 280.460 + 0.9856474 * n  # 평균 경도 L (0~360도)
    while L > 360:  # 경도가 360도가 넘어갈때 보정하는 식
        L = L - 360

    g = 357.528 + 0.9856003 * n  # 평균 근점이각 g (0~360도)
    while g > 360:  # 평균 근점이각이 360도가 넘어갈때 보정하는 식
        g = g - 360

    I = L + 1.915 * math.sin(math.radians(g)) + 0.02 * math.sin(2 * math.radians(g))  # 황도 경도 I(대문자 i) (0~360도)
    while I > 360:  # 황도 경도가 360도가 넘어갈때 보정하는 식
        I = I - 360

    ep = 23.439 - 0.0000004 * n

    ra = 360 + math.atan2(math.cos(math.radians(ep)) * math.sin(math.radians(I)),
                          math.cos(math.radians(I))) * 180 / np.pi  # 적경 ra (0~360도)
    while ra > 360:  # 적경이 360도가 넘어갈때 보정하는 식
        ra = ra - 360

    dec = math.asin(math.sin(math.radians(ep)) * math.sin(math.radians(I))) * 180 / np.pi  # 적위 dec

    gmst = 6.697375 + 0.0657098242 * n + hour  # Greenwich 평균항성시 (각도)
    gmst = gmst * 15
    while gmst > 360:  # Greenwich 평균항성시를 24시간으로 표현
        gmst = gmst - 360

    gmst = gmst / 15

    lmst = gmst + lte / 15  # 레일 계측 지점의 평균항성시간 / lte(longitude) : 경도 / (예 : 서울지방 경도 : 126.98333)
    ha = lmst - ra / 15  # 시간각 (-12<ha<12)
    ha = ha * 15  # 시간각(24시간)을 각도로 표현

    el = math.asin(math.sin(math.radians(dec)) * math.sin(math.radians(lat)) + math.cos(math.radians(dec)) * math.cos(
        math.radians(lat)) * math.cos(math.radians(ha))) * 180 / np.pi  # 태양고도각 el / lat(latitude) : 위도
    az = math.asin(-math.cos(math.radians(dec)) * math.sin(math.radians(ha)) / math.cos(
        math.radians(el))) * 180 / np.pi  # 태양 방위각 az

    if math.sin(math.radians(el)) >= math.sin(math.radians(dec)) / math.sin(math.radians(lat)):
        az = 180 - az;

    elif math.sin(math.radians(el)) <= math.sin(math.radians(dec)) / math.sin(math.radians(lat)) and ha > 0:
        az = 360 + az;

    if az < 0:  # 태양 방위각 보정식
        az = az + 360
    else:
        while az > 360:
            az = az - 360

    return az,el

def TSI(time_raw):
    per_date = datetime.datetime(time_raw.year, 1, 2)
    seq_date = time_raw
    day_delta = seq_date - per_date
    day_delta_raw = day_delta.days
    theta = (day_delta_raw * (360 / orbit_days)) * (np.pi / 180)
    r = (au * (1 - pow(e, 2))) / (1 + e * math.cos(theta)) / (pow(10, 8))
    TSI_raw = solar_c * pow((1.4957911194950864 / r), 2)
    return TSI_raw
