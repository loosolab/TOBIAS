#!/usr/bin/env python3

from argparse import ArgumentParser
import logging
import boto3
from botocore.client import Config
import botocore
import yaml
import os
import sys
import fnmatch

logging.basicConfig(format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s', level=logging.INFO)

#--------------------------------------------------------------------------------------------------------#

def read_config_yaml(config_yaml):
    
    #read yaml file
    with open(config_yaml, 'r') as stream:
        try:
            yaml_input = yaml.safe_load(stream)
        except yaml.YAMLError as exc:
            logging.error(exc)
            sys.exit(1)

    #Check required Yaml keys
    if not all(key in yaml_input for key in ['endpoint', 'buckets']):
        logging.error('Missing required_keys \'buckets\' and  \'endpoint\'')
        sys.exit(1)

    if not all(key in yaml_input for key in ['username', 'accesskey']):
        yaml_input['username'] = ''
        yaml_input['accesskey'] = ''
    return yaml_input

#--------------------------------------------------------------------------------------------------------#

def s3_downloader(client, bucket, download_list):
    for key in download_list:
        #Create folder Path if dosen't exists
        if not os.path.exists('S3_download/' + bucket + '/' + os.path.dirname(key)):
            logging.info('Create directory: ' + bucket + '/' +  os.path.dirname(key))
            os.makedirs('S3_download/' + bucket + '/' + os.path.dirname(key))
        #Check if file already exists
        if os.path.isfile('S3_download/' + bucket + '/' + key):
            logging.warning('File: S3_download/' + bucket + '/' + os.path.basename(key) + ' already exists!')
            continue
        #Download file
        client.download_file(bucket, key, 'S3_download/' + bucket + '/' + key)

#--------------------------------------------------------------------------------------------------------#

def s3_client(yaml_input):
    
    bsession = boto3.Session()
    #Create client 
    client = bsession.client('s3',
            endpoint_url = yaml_input['endpoint'], 
            aws_access_key_id = yaml_input['username'],
            aws_secret_access_key = yaml_input['accesskey'])

    #read filter files for download of each 
    for bucket in yaml_input['buckets']:
        #try to access bucket with if credentials if no credentials set boto search for credentails under ~/.aws/credentials
        try:
            #if no credentials set skip login and go to public access
            if len(yaml_input['username']) == 0:
                raise UnboundLocalError()
            bucket_objects = [obj['Key'] for obj in client.list_objects_v2(Bucket = bucket)['Contents']]
            logging.info('Connect to bucket: ' + str(bucket))
        except Exception as e:
            if len(yaml_input['username']) > 0:
                logging.warning('Cloud not access bucket ' + str(bucket) + ' with credentials')
            #if credentails dosent work try to access bucket as public bucket            
            try:
                client.meta.events.register('choose-signer.s3.*', botocore.handlers.disable_signing)#disable sigin with credentials for public buckets
                bucket_objects = [obj['Key'] for obj in client.list_objects_v2(Bucket = bucket)['Contents']]
                logging.info('Connect to public bucket: ' + str(bucket))
            except Exception as e:
                logging.error('Cloud not access bucket: ' + str(bucket) + ' ' + str(e))
                continue
            
        #Download files 
        #no pattern downlaod whole Bucket
        if not yaml_input['buckets'][bucket] or len(yaml_input['buckets'][bucket]) == 0:
            s3_downloader(client, bucket, bucket_objects)
        else:
            for pattern in yaml_input['buckets'][bucket]:
                #Check if bucket has files which match the pattern
                if len(fnmatch.filter(bucket_objects, pattern)) == 0:
                    logging.warning('Cloud not find file for pattern: ' + pattern)
                    continue
                #download files which machtes the pattern
                s3_downloader(client, bucket, fnmatch.filter(bucket_objects, pattern))
                

#--------------------------------------------------------------------------------------------------------#

def argparsser():
    parser = ArgumentParser()
    parser.add_argument("--yaml_file", dest="yamfile", type=str,   help="Config Yaml file.")
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    args = parser.parse_args()
    return args

#--------------------------------------------------------------------------------------------------------#

def main():
    args = argparsser()
    s3_client(read_config_yaml(args.yamfile))

#--------------------------------------------------------------------------------------------------------#

if __name__ == '__main__':
    main()