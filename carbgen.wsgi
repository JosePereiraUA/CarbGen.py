import os, sys

PROJECT_DIR = '/home/ssantos/venv/carbgen'
sys.path.insert(0, PROJECT_DIR)


def execfile(filename):
    globals = dict( __file__ = filename )
    exec( open(filename).read(), globals )

activate_this = os.path.join( PROJECT_DIR, 'venv/bin', 'activate_this.py' )
execfile( activate_this )


from app import app as application
from app import db
db.create_all()
