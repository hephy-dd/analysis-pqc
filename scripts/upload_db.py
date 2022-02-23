import subprocess
import logging
import pdb
import os
import argparse
import glob

logger = logging.getLogger(__name__)

logging.basicConfig(filename='xml_upload.log', level=logging.INFO)#, format='%(levelname)s:%(name)s:%(message)s')

script_path = os.path.dirname(os.path.realpath(__file__))  # scripts directory


def get_run_number(db_instance):
    ext_filename = os.path.join(script_path, 'ext', 'rhapi.py')    
    p1 = subprocess.run( ['python3', ext_filename, '--login', '--url=https://cmsdca.cern.ch/trk_rhapi', '-f', 'csv', '--clean',f"select r.run_number from trker_{db_instance}.trk_ot_test_nextrun_v r"], capture_output=True)
    #python3 ext/rhapi.py --login --url=https://cmsdca.cern.ch/trk_rhapi -f csv --clean "select r.run_number from trker_cmsr.trk_ot_test_nextrun_v r"
    answer = p1.stdout.decode()
    answer = answer.split()

    #return only the number (run_number) which corresponds to the output message
    logging.info(f'Query run-number replied with {answer}')
    return answer[1]


def upload_to_db(filename, db_instance):
    ext_filename = os.path.join(script_path, 'ext', 'cmsdbldr_client.py')    
    try:
        p1 = subprocess.run(['python3', ext_filename, '--login', f'--url=https://cmsdca.cern.ch/trk_loader/trker/{db_instance}', f'{filename}'], capture_output=True)
        # python3 ext/cmsdbldr_client.py --login --url=https://cmsdca.cern.ch/trk_loader/trker/int2r filename
        answer = p1.stdout.decode()
        answer = answer.split()
        return answer
    except Exception as error:
        print(error)


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

    
def run(path,db):
   
    if db == 'development': db_instance = 'int2r'
    elif db == 'production': db_instance = 'cmsr'
    else: raise Exception('DB not understood, choose production or development')

    filenames=[]
    for root, dirs, files in os.walk(path):
        xml_files=[filename for filename in files if os.path.splitext(filename)[1] == '.xml']
        for xml_file in xml_files: filenames.append(os.path.join(root,xml_file))
    number_of_files=len(filenames)
    print(f'Upload {number_of_files} files to {db} database, continue? (yes)')
    choice = input().lower()

    uploaded=0
    failed=[]
    if choice == 'yes':
        logging.info('-----------------------------------------------------------------')
        logging.info(f'Uploading {number_of_files} files from {path} to {db} database\n')
        for ifile,filename in enumerate(filenames):

            ##Query of run number not necessary, will be assigned on upload
            #run=get_run_number(db_instance) 
            #update_file(filename,'<RUN_NUMBER></RUN_NUMBER>',f'<RUN_NUMBER>{run}</RUN_NUMBER>')
            print('Uploading file',ifile,end='\r')
            answer=upload_to_db(filename,db_instance)
            if answer[3]!='200':
                logging.info(f'Failed: Uploading file {filename} to {db} database, replied with {answer}')
                failed.append(filename)
            else:
                uploaded+=1
        print('\n')
        print(uploaded,'files uploaded')
        logging.info(f'Uploaded/Total: {uploaded}/{number_of_files}\n')

    else:
        print('Nothing uploaded')    

if __name__ == "__main__":

    args = parse_args()
    path = args.path
    db = args.database

    run(path,db)
