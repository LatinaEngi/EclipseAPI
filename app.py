#!flask/bin/python
from flask import Flask, jsonify, make_response, request, abort

from astropy.time import Time
from astropy.coordinates import EarthLocation, get_sun, AltAz

from moon import moon_distance, get_umbra, estimate_umbra

app = Flask(__name__)

@app.route('/', methods=['GET'])
def say_hello():
    return jsonify(message='Hello! Welcome to the 2017 Eclipse App.')

@app.route('/results', methods=['POST'])
def do_calculations():
    if not request.json or not 'location' in request.json or not 'startTime' in request.json or not 'endTime' in request.json:
        abort(400)

    app_vars = {
                'location' : 'Nashville, TN',
                'c2time' : '18:27:03.9',
                'c3time' : '18:29:43.0',
                'd_est' : '',
                'd_real' : '',
                'su_est' : '',
                'su_real' : '',
                'lat' : '',
                'long' : '',
                'sun_alt' : '',
                'percent_diff' : ''
    }

    app_vars['location'] = request.json['location']
    app_vars['c2time'] = request.json['startTime']
    app_vars['c3time'] = request.json['endTime']

    loc = EarthLocation.of_address(app_vars['location'])
    t2 = Time('2017-08-21 ' + app_vars['c2time'])
    t3 = Time('2017-08-21 ' + app_vars['c3time'])

    dm, reald = moon_distance(loc, t2, t3)

    sun = get_sun(t2)
    sun_altaz = sun.transform_to(AltAz(obstime=t2, location=loc))

    su_est = estimate_umbra(loc, t2, t3)
    su_real = get_umbra(loc, t2)

    app_vars['d_est'] = '{:.7}'.format(dm)
    app_vars['d_real'] = '{:.7}'.format(reald)
    app_vars['su_est'] = '{:.4}'.format(su_est)
    app_vars['su_real'] = '{:.4}'.format(su_real)
    app_vars['lat'] = '{:.4}'.format(loc.latitude.value)
    app_vars['long'] = '{:.4}'.format(loc.longitude.value)
    app_vars['sun_alt'] = '{:.4}'.format(sun_altaz.alt.value)
    app_vars['percent_diff'] = '{:.2}'.format((dm - reald) / reald * 100.)

    return jsonify(app_vars)

@app.errorhandler(404)
def not_found(error):
    return make_response(jsonify({'error': 'Not found'}), 404)

if __name__ == '__main__':
    app.run(debug=False)
