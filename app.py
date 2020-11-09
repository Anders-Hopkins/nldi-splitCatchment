## StreamStats delineation script flask app server wrapper

# -----------------------------------------------------
# Martyn Smith USGS
# 09/30/2019
# StreamStats Delineation script flask server
#
# run with: "python -m flask run"
#
# -----------------------------------------------------

from flask import Flask, request, jsonify
from flask_cors import CORS, cross_origin
from datetime import datetime
import delineate

app = Flask(__name__)
cors = CORS(app)
app.config['CORS_HEADERS'] = 'Content-Type'

@app.route("/")

def home():
    return "sample delineation server"

@app.route("/delineate")
@cross_origin(origin='*')
def main():

    region = request.args.get('region')
    lat = float(request.args.get('lat'))
    lng = float(request.args.get('lng'))

    print(region,lat,lng)
    dataPath = 'C:/NYBackup/GitHub/ss-delineate/data/'

    #start main program
    results = delineate.Watershed(lat,lng,region,dataPath)
    return jsonify(results.serialize())