from flask_sqlalchemy import SQLAlchemy
from flask_wtf.csrf import CSRFProtect
from flask import Flask


from config import Config


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Flask APP Instantiation
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
app = Flask(__name__)
app.config.from_object(Config)

csrf = CSRFProtect(app)
db = SQLAlchemy(app)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# Session Scope context manager
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
from sqlalchemy.orm import sessionmaker
from sqlalchemy import create_engine
from contextlib import contextmanager

engine = create_engine(Config.SQLALCHEMY_DATABASE_URI)
Session = sessionmaker(bind=engine)

@contextmanager
def session_scope():
    '''
    Provide a transactional scope around a series of operations.

    Ref. http://docs.sqlalchemy.org/en/latest/orm/session_basics.html
    '''
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.close()



from app import routes


if __name__ == '__main__':
    db.create_all()
    app.run()
    # app.run(host='0.0.0.0', port=5000, debug=True)
