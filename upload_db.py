import subprocess
import logging
import pdb
import os
import argparse
import glob

logger = logging.getLogger(__name__)

logging.basicConfig(filename='xml_upload.log', level=logging.INFO)#, format='%(levelname)s:%(name)s:%(message)s')


def get_run_number():
    
    p1 = subprocess.run( ['python3', 'rhapi.py', '--login', '--url=https://cmsomsdet.cern.ch/tracker-resthub', '-f', 'csv', '--clean',"select r.run_number from trker_cmsr.trk_ot_test_nextrun_v r"], capture_output=True)
    
    answer = p1.stdout.decode()
    answer = answer.split()

    #return only the number (run_number) which corresponds to the output message
    logging.info(f'Query run-number replied with {answer}')
    return answer[1]


def upload_to_db(filename, db_instance):

    try:
        p1 = subprocess.run(['python3', 'cmsdbldr_client.py', '--login', f'--url=https://cmsdca.cern.ch/trk_loader/trker/{db_instance}', '{filename}'],  capture_output=True)
        ###p1 = subprocess.run(
             #   'python cmsdbldr_client.py --login --url=https://cmsdca.cern.ch/trk_loader/trker/cmsr {}',format(file),
              #  capture_output=True)

        answer = p1.stdout.decode()
        answer = answer.split()
        logging.info(f'Uploading file {filename} to {db} database, replied with {answer}')
    except Exception as error:
        print(error)
        #logging.exception('Exception occured: {}'.format(error))


def parse_args():

    parser = argparse.ArgumentParser()
    parser.add_argument('path')
    parser.add_argument('-db', dest='database', default='production', help='choose DB: production or development')

    return parser.parse_args()


def update_file(filename,pattern,replace):
    with open(filename,'r') as file_in:
        lines=file_in.read()
        lines=lines.replace(pattern,replace)
    with open(filename,'w') as file_out:
        file_out.write(lines)
    logging.info(f'Updating file {filename}')


def find_files_ext(filenames,ext):
    '''
    filenames: str, list of files
    ext: str, e.g. '.xml'
    return: reducced list of files
    '''
    filenames_out=[]
    for filename in filenames:
        _, extension = os.path.splitext(filename)
        is_xml_file = extension == ext
        if is_xml_file:
            filenames_out.append(filename)
    return filenames_out

    
def run(path,db):

    if db == 'development': db_instance = 'int2r'
    elif db == 'production': db_instance = 'cmsr'
    else: raise Exception('DB not understood, choose production or development')

    filenames=[]
    
    for root, dirs, files in os.walk(path):
        xml_files=find_files_ext(files,'.xml')
        for xml_file in xml_files: filenames.append(root + os.sep + xml_file)

    print(f'Upload {len(filenames)} files to {db} database, continue? (yes)')
    choice = input().lower()
    if choice == 'yes':
        for ifile,filename in enumerate(filenames):
            if ifile>0: break
        
            #run=get_run_number() ##Query of run number not necessary, will be assigned on upload
            #update_file(filename,'<RUN_NUMBER></RUN_NUMBER>',f'<RUN_NUMBER>{run}</RUN_NUMBER>')

            upload_to_db(filename,db_instance)
    else:
        print('Nothing uploaded')    

if __name__ == "__main__":

    args = parse_args()
    path = args.path
    db = args.database

    run(path,db)
