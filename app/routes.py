from flask import request, render_template, url_for
from flask import Flask, Response

from app import app, db, session_scope, Config
from models import Job, Statistics

import urllib2
import time
import json
import os

import carbgen
#--------------------------------------------------------------------
# HELPER FUNCTIONS - HELPER FUNCTIONS - HELPER FUNCTIONS - HELPER FUN
#--------------------------------------------------------------------
def run_job(job_id, parameters):
    
    # print url_for('static', filename='treta.ff')
    ff = '/static/sheet_ga.ff'
    sheet_gaff_dir = Config.BASE_DIR + os.path.abspath(ff)
    
    zip_content = carbgen.batch_create(
        int(parameters['batch_size']),
        "residue",              # the name of the residue!
        int(parameters['x']),
        int(parameters['y']),
        int(parameters['z']),
        0.8 - 0.02 * float(parameters['pf']),
        int(parameters['co']),
        int(parameters['coo']),
        int(parameters['ror']),
        int(parameters['ch']),
        sheet_gaff_dir
    )

    # emulate a long runnning process
    time.sleep(2)
    with session_scope() as session:
        job = session.query(Job).get(job_id)
        if job is not None:
            # and update its status
            job.result = buffer(zip_content)
            job.completed = True
            session.commit()


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@app.before_first_request
def initialize():
    '''
    Provides the application context with a pool
    of workers to build residues on the background.
    '''
    from multiprocessing.pool import ThreadPool

    app.logger.info('Creating a ThreadPool obect with %d workers', Config.MAX_POOL_SIZE)
    app.pool = ThreadPool(Config.MAX_POOL_SIZE)



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @app.teardown_appcontext
# def close_connection(exception):
#     '''
#     Utility callback that closes the DB connection on app exit.
#     '''
#     db = getattr(g, '_database', None)
#     if db is not None:
#         db.close()
#     app.logger.info('Closed db connection')



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@app.route('/')
def input():
    '''
    Main entry point to the application.
    The template is rendered and a CSRF token
    is appended to prevent cross-site request forgery attacks.
    '''
    
    return render_template('index.html')




#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@app.route('/generate', methods=['POST'])
def generate():
    
    input_data = {k:v for k,v in request.form.iteritems()}
    app.logger.info('Received new job request with parameters: %s', input_data)
    print input_data
    if not carbgen.validate_parameters(**input_data):
        # status = 422    # Unprocessable Entity
        # msg = '{"msg": "Unable to fulfill request"}'
        return render_template('index.html')

    try:
        # generate a new job
        job = Job()
        db.session.add(job)
        db.session.commit()

        # and assign it to the worker pool
        # run_job(job.id, input_data)
        app.pool.apply_async(run_job, args=(job.id, input_data))
    
    except Exception as e:
        app.logger.error('Exception %s', e.message)
        msg = '{"msg": "Unable to fulfill request"}'
        status = 500

    else:
        app.logger.info('Created new job %s', job)
        msg = '{"job": %d}' % job.id
        status = 200
    
    return Response(msg, status=status, mimetype='application/json')

    # # Create new JOB
    # job = Job()
    # db.session.add(job)
    # db.session.commit()
    # app.logger.info('Created new job %s', job)

    # try:
    #     # assign the job to the worker pool
    #     # generate(job_id, input_data)

    #     app.pool.apply_async(generate, args=(job.id, input_data))
    
    # except Exception as e:
    #     app.logger.error('%s', e.message)
    #     status = 422    # Unprocessable Entity
    #     msg = '{"msg": "Unable to fulfill request"}'
    # else:
    #     status = 200
    #     msg = '{"job": %d}' % job.id
    
    # # return a message to the client
    # return Response(msg, status=status, mimetype='application/json')



#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@app.route('/pool_result', methods=['POST'])
def pool_result():
    '''
    Endpoint to check if a job with a given ID (in a POST request)
    is ready for download. 
    '''

    # get the job ID from the POST data
    job_id = int(request.form.get('job'))
    
    job = Job.query.get(job_id)

    # validate the query output
    if job is None:
        # if the job is None, then the jod ID is unknown (ERROR)
        app.logger.error('Unknown job %s', job_id)
        msg = '{"msg": "Unknown job id"}'
        status = 400

    else:
        # otherwise, get the completion status from the query output
        app.logger.debug('status for job %s', job)
        msg = '{"completed": %s}' % str(job.completed).lower()
        status = 200
    
    # return the Response object
    return Response(msg, status=status, mimetype='application/json')
    


#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
@app.route('/download_result', methods=['GET'])
def download_result():
    job_id = int(request.args.get('job'))
    
    app.logger.debug('Downloading result for job %s', job_id)
    
    job = Job.query.get(job_id)
    if job is None:
        app.logger.error('Unknown job %s', job_id)
        response = Response(
            '{"msg": "Unknown job id"}',
            mimetype='application/json',
            status=400)

    else:
        # get result from row tuple (zipped-file content)
        # the result is a python's buffer object
        # that needs to be converted to a str object 
        result = str(job.result)
        

        # delete entry from database
        db.session.delete(job)
        #db.session.commit()
        app.logger.debug('Deleted job %s', job)

        try:
            client_ip = request.environ.get('HTTP_X_REAL_IP', request.remote_addr)
            url = 'http://ipinfo.io/%s/json' % client_ip
            response = urllib2.urlopen(url)
            client_data = json.load(response)

            
        except Exception as e:
            app.logger.fatal(e, exc_info=True)
        else:
            country = client_data.get('country', 'UNK')
            stats = Statistics.query.get(country)
            if stats is None:
                stats = Statistics(country)
                db.session.add(stats)
            stats.count += 1
        db.session.commit()



        # construct response
        response = Response(
            result,
            mimetype='application/zip',
            headers={
                'Content-Disposition': 'attachment;filename=carbgen.zip'
            })
    
    return response








#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# LEGACY ENDPOINTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# @app.route('/confirm', methods=['POST'])
# def confirm():
    
#     input_data = fixDictionary(request.form.to_dict(flat=False))
#     input_data['pf'] = 0.8 - input_data['pf'] * 0.02
    
#     zip_content = residue.batch_create(
#         input_data['batch_size'],
#         input_data['name'],
#         input_data['x'],
#         input_data['y'],
#         input_data['z'],
#         input_data['pf'],
#         input_data['co'],
#         input_data['coo'],
#         input_data['ror'],
#         input_data['ch']
#     )
    
#     # return zipped files
#     zipname = input_data['name'] + 'zip'
#     return Response(zip_content,
#             mimetype='application/zip',
#             headers={
#                 'Content-Disposition': 'attachment;filename=%s'%zipname
#             }
#         )


# def fixDictionary(dictionary):
#     for i in dictionary:
#         if not i == "name": dictionary[i] = int(dictionary[i][0])
#         else: dictionary[i] = str(dictionary[i][0])
#     return dictionary



# if __name__ == '__main__':
#     # CREATE TABLE JOBS (
#     #    ID           INTEGER PRIMARY KEY AUTOINCREMENT,
#     #    COMPLETED    INTEGER NOT NULL,
#     #    RESULT       BLOB
#     # );
#     # http://flask.pocoo.org/docs/0.12/patterns/sqlite3/

#     # with app.app_context():
#     #     db = get_db()
#     #     with app.open_resource('schema.sql', mode='r') as f:
#     #         db.cursor().executescript(f.read())
#     #     db.commit()
#     app.run(host='0.0.0.0')
    
