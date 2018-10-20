import re
from utils import *
import itertools
import csv
import sys
import ConfigParser
import os
from io import StringIO
from passlib.hash import argon2
from bson.json_util import dumps
from json import dumps
 
#flask import
from flask import Flask
from flask import session
from flask.ext.session import Session
from flask import Response
from flask import stream_with_context
from flask import request
from flask import make_response
from flask import request
from flask import send_file
from flask import g
from flask import redirect
from flask import url_for
from flask import abort
from flask import render_template
from flask import flash
from flask import jsonify
from flask import send_from_directory
from flask.ext.compress import Compress
from flask.ext.runner import Runner
from flask_errormail import mail_on_500
from flask_debugtoolbar import DebugToolbarExtension 
from flask.ext.cache import Cache
from flask import Flask, jsonify, request
from flask_jwt_extended import ( JWTManager, jwt_required, create_access_token, get_jwt_identity)
import logging
logging.getLogger().addHandler(logging.StreamHandler())
logging.getLogger().setLevel(logging.INFO)
# define app
# Load default config and override config from an environment variable
from config import config
if config.LOCAL:
    print 'LOCAL'
    app = Flask(__name__,static_url_path='/static', static_folder='../static', template_folder='../templates')
    app.config.from_pyfile('../local.cfg')
else:
    print 'SERVER'
    app = Flask(__name__, template_folder='../templates')
    app.config.from_pyfile('../phenopolis.cfg')
app.config['JWT_SECRET_KEY'] = 'super-secret'  # Change this!
jwt = JWTManager(app)
Compress(app)
cache = Cache(app,config={'CACHE_TYPE': 'simple'})
from functools import wraps 
def requires_auth(f):
    @wraps(f)
    def decorated(*args, **kwargs):
        if session:
          if 'user' in session: 
             return f(*args, **kwargs)
        if request.method == 'POST':
          username=request.form['user']
          password=request.form['password']
          if check_auth(username,password):
             return f(*args, **kwargs)
        print 'Not Logged In - Redirect to home to login'
        if config.LOCAL:
           return redirect('/')
    return decorated
# neo4j import
from neo4j.v1 import GraphDatabase, basic_auth, CypherError
import neo4j_setup

uri='bolt://0.0.0.0:7687'
password=1
neo4j_driver = GraphDatabase.driver(uri, auth=basic_auth("neo4j", password))
neo4j_session=neo4j_driver.session()

# views
import views.endpoints
# work in progress, comment out if not needed



