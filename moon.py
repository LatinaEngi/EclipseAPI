# Function to calculate the distance to the moon

import requests
import pandas as pd
import numpy as np
from astropy.constants import R_earth, R_sun
from astropy import units as u
from astropy.coordinates import get_sun, get_moon, EarthLocation, AltAz
from bs4 import BeautifulSoup


def moon_distance(loc, t2, t3):
    # Get the distance to the moon via the umbra diameter estimate

    if isinstance(loc, type(u'')) or isinstance(loc, type('')):
        loc = EarthLocation.of_address(loc)

    # Location of Sun & Moon
    sun = get_sun(t2)
    moon = get_moon(t2, loc)

    # Umbra size, from time and location
    su = estimate_umbra(loc, t2, t3)

    # Constants in km
    re = R_earth.to(u.km).value
    ds = sun.distance.to(u.km).value
    rs = R_sun.to(u.km).value
    rm = 1738.1

    dm = (2 * rm - su) * ds / (2 * (rs - rm)) + re
    #print('Distance estimate: {:.7} km (Actual: {:.7}). Ratio: {:.2}. Difference: {:.7}'.format(dm, moon.distance, moon.distance / (dm * u.km), dm * u.km - moon.distance))

    # moon.distance is from surface, not center
    return dm, moon.distance.to(u.km).value + re


def estimate_umbra(loc, t2, t3, verbose=False):
    # Estimate the size of the umbra from the contact times (not 100% sure on this yet)
    delta_t = t3 - t2

    # Location of Sun
    sun = get_sun(t2)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    # Umbra size, from time and location
    s_umbra = 2 * np.pi * R_earth * np.cos(loc.latitude.to(u.rad)) / u.day * delta_t
    s_umbra /= np.tan(sun_altaz.alt.to(u.rad))  # correction for altitude (correct?)
    su = s_umbra.to(u.km).value

    # Extra corrections for projection/location from central line
    # Need to check if this is the right approach
    df = get_path()
    df['sep'] = ang_sep(loc, df['clon'], df['clat']) * 2 * np.pi * R_earth.to(u.km).value / 360.
    ind = df[df['sep'] == min(df['sep'])].index[0]

    # Replace with finer version
    x2 = 41.
    x1 = 0.
    df_fine = pd.DataFrame(np.arange(x1, x2, 1), columns=['xval'])
    df_fine['clat'] = make_fine('clat', ind, df, x1=x1, x2=x2)
    df_fine['clon'] = make_fine('clon', ind, df, x1=x1, x2=x2)
    df_fine['width'] = make_fine('width', ind, df, x1=x1, x2=x2)

    # Get new best index
    df_fine['sep'] = ang_sep(loc, df_fine['clon'], df_fine['clat']) * 2 * np.pi * R_earth.to(u.km).value / 360.
    ind = df_fine[df_fine['sep'] == min(df_fine['sep'])].index[0]
    df = df_fine

    # Get the angle
    factor = np.pi / 180
    lat1 = df.loc[ind, 'clat'] * factor
    lat2 = lat1
    dlon = (df.loc[ind, 'clon'] - df.loc[ind - 1, 'clon']) * factor
    theta1 = np.arctan2(np.sin(dlon) * np.cos(lat2),
                        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)) / factor
    lat2 = df.loc[ind - 1, 'clat'] * factor
    theta2 = np.arctan2(np.sin(dlon) * np.cos(lat2),
                        np.cos(lat1) * np.sin(lat2) - np.sin(lat1) * np.cos(lat2) * np.cos(dlon)) / factor
    theta = abs(theta1 - theta2)

    # Modify the umbra size
    # This still needs more work as it doesn't give an approximately correct value
    ratio = 0.5*df.loc[ind, 'width'] / (0.5*df.loc[ind, 'width']-df.loc[ind, 'sep'])
    su2 = su * ratio / np.cos(theta * factor)

    return su2


def get_umbra(loc, t2):
    # Get the real size of the umbra

    sun = get_sun(t2)
    moon = get_moon(t2, loc)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    # Constants in km
    re = R_earth.to(u.km).value
    ds = sun.distance.to(u.km).value
    rs = R_sun.to(u.km).value
    rm = 1738.1

    alpha = np.arcsin((rs - rm) / (ds - moon.distance.to(u.km).value))
    su = 2 * (rm / np.sin(alpha) - (moon.distance.to(u.km).value - re)) * np.tan(alpha)
    su / np.tan(sun_altaz.alt.to(u.rad))  # not 100% sure of this

    return su


def convert_to_decimal(e1, e2):
    # Quick function to convert the inputs to float
    if e2[-1] == 'N' or e2[-1] == 'E':
        sign = 1.
    else:
        sign = -1.

    return (float(e1) + float(e2[:-1]) / 60) * sign

def get_path():
    # Get the path of the eclipse from the NASA page
    
    # url = 'https://eclipse.gsfc.nasa.gov/SEpath/SEpath2001/SE2017Aug21Tpath.html'
    #r = requests.get(url)
    #html = r.content
    #r.close()
    
    r = open('SE2017Aug21Tpath.html', 'r', encoding='latin-1')
    html = r.read()
    r.close()

    # Extract the data
    soup = BeautifulSoup(html, 'html.parser')
    raw_data = soup.find('pre').get_text()

    parse_data = list()
    for i, raw_line in enumerate(raw_data.split('\r')):
        line = raw_line.strip()
        if line == '' or len(line) == 0 or line.startswith('Limits') or i < 7:
            continue

        elems = line.split()
        if len(elems) < 17:
            continue

        utc = elems[0]
        nlat = convert_to_decimal(elems[1], elems[2])
        nlon = convert_to_decimal(elems[3], elems[4])
        slat = convert_to_decimal(elems[5], elems[6])
        slon = convert_to_decimal(elems[7], elems[8])
        clat = convert_to_decimal(elems[9], elems[10])
        clon = convert_to_decimal(elems[11], elems[12])
        ms_ratio = float(elems[13])
        sun_alt = float(elems[14])
        sun_az = float(elems[15])
        width = float(elems[16])
        dur = elems[17]

        parse_data.append([utc, nlat, nlon, slat, slon, clat, clon, ms_ratio, sun_alt, sun_az, width, dur])

    df = pd.DataFrame(parse_data,columns=['time', 'nlat', 'nlon', 'slat', 'slon', 'clat', 'clon', 'msratio', 'alt', 'az', 'width', 'dur'])

    return df

def ang_sep(loc, ra1, dec1):
    # Using Vicenty Formula (http://en.wikipedia.org/wiki/Great-circle_distance)
    # and adapting from astropy's SkyCoord
    long = loc.longitude.value
    lat = loc.latitude.value
    factor = np.pi / 180
    sdlon = np.sin((long - ra1) * factor)  # RA is longitude
    cdlon = np.cos((long - ra1) * factor)
    slat1 = np.sin(dec1 * factor)  # Dec is latitude
    slat2 = np.sin(lat * factor)
    clat1 = np.cos(dec1 * factor)
    clat2 = np.cos(lat * factor)

    num1 = clat2 * sdlon
    num2 = clat1 * slat2 - slat1 * clat2 * cdlon
    numerator = np.sqrt(num1 ** 2 + num2 ** 2)
    denominator = slat1 * slat2 + clat1 * clat2 * cdlon

    # Output in degrees
    return np.arctan2(numerator, denominator) / factor


def make_fine(col, ind, df, x1=0, x2=21):
    # Create a finer version of the path
    y2 = df.loc[ind+1, col]
    y1 = df.loc[ind-1, col]
    m = (y2-y1)/(x2-x1)
    return np.arange(x1, x2, 1) * m + y1
