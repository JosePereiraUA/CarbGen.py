import os

basedir = os.path.abspath(os.path.dirname(__file__))

dbdir = os.path.join(os.path.dirname(basedir), 'persistence')

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-guess'
    SQLALCHEMY_DATABASE_URI = os.environ.get('DATABASE_URL') or \
        'sqlite:///' + os.path.join(dbdir, 'carbgen.db')
    SQLALCHEMY_TRACK_MODIFICATIONS = False
    MAX_POOL_SIZE = 2
    BASE_DIR = basedir
    DB_DIR = dbdir
