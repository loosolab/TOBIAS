#!/usr/bin/env python

"""
DownloadData: Download test data from the loosolab s3 server

@authors: Mette Bentsen and Philipp Goymann
@contact: mette.bentsen (at) mpi-bn.mpg.de
@license: MIT
"""

import boto3
from botocore.client import Config
import botocore
import yaml
import os
import sys
import fnmatch

from tobias.utils.logger import TobiasLogger
from tobias.utils.utilities import make_directory

#--------------------------------------------------------------------------------------------------------#
def read_config_yaml(config_yaml, logger):
	
	#Read yaml file
	with open(config_yaml, 'r') as stream:
		try:
			config_dict = yaml.safe_load(stream)
		except yaml.YAMLError as exc:
			logger.error(exc)
			sys.exit(1)

	#Check required yaml keys
	if not all(key in config_dict for key in ['endpoint', 'buckets']):
		logger.error('Missing required_keys \'buckets\' and  \'endpoint\'')
		sys.exit(1)

	if not all(key in config_dict for key in ['username', 'accesskey']):
		config_dict['username'] = None
		config_dict['accesskey'] = None

	return config_dict

#--------------------------------------------------------------------------------------------------------#
def s3_downloader(client, bucket, download_list, logger, force):
	""" Target is the name of the folder to save files to """

	for s3_file in download_list:

		#Create folder Path if doesn't exists
		target_path = os.path.join(bucket, s3_file) #if target == None else os.path.join(target, s3_file)

		#Create directory if it is not already there
		if not os.path.exists(os.path.dirname(target_path)):
			logger.info('Creating directory: {0}'.format(os.path.dirname(target_path)))
			make_directory(os.path.dirname(target_path))

		#Check if file already exists
		if os.path.isfile(target_path):
			if force == False:	#if file exists and force is off
				logger.warning("File '{0}' already exists! Use --force to overwrite.".format(target_path))
				continue

		#Download file
		logger.info("Downloading file '{0}' to target '{1}'".format(s3_file, target_path))
		client.download_file(bucket, s3_file, target_path)

#--------------------------------------------------------------------------------------------------------#
def s3_client(config_dict, logger, force):
	
	bsession = boto3.Session()

	#Create client 
	client = bsession.client('s3',
			endpoint_url = config_dict['endpoint'], 
			aws_access_key_id = config_dict['username'],
			aws_secret_access_key = config_dict['accesskey'])

	#Read filter files for download of each 
	for bucket in config_dict['buckets']:

		#Try to access bucket with if credentials if no credentials set boto search for credentails under ~/.aws/credentials
		try:
			#if no credentials set skip login and go to public access
			if config_dict['username'] is None:
				raise UnboundLocalError()

			bucket_objects = [obj['Key'] for obj in client.list_objects_v2(Bucket = bucket)['Contents']]
			logger.info('Connected to bucket: ' + str(bucket))

		except Exception as e:
			
			#If credentials were set, but somehow did not connect
			if config_dict['username'] != None:
				logger.warning('Could not not access bucket ' + str(bucket) + ' with credentials')

			#If credentails doesn't work, try to access bucket as public bucket            
			try:
				client.meta.events.register('choose-signer.s3.*', botocore.handlers.disable_signing)#disable sigin with credentials for public buckets
				bucket_objects = [obj['Key'] for obj in client.list_objects_v2(Bucket = bucket)['Contents']]
				logger.info('Connected to public bucket: ' + str(bucket))

			except Exception as e:
				logger.error('Could not access bucket: ' + str(bucket) + ' ' + str(e))
				continue
			
		#Download files 
		#If there was no pattern; download whole Bucket
		if not config_dict['buckets'][bucket] or len(config_dict['buckets'][bucket]) == 0:
			s3_downloader(client, bucket, bucket_objects)
		else:
			for pattern in config_dict['buckets'][bucket]:

				#Check if bucket has files which match the pattern
				pattern_files = fnmatch.filter(bucket_objects, pattern)
				if len(pattern_files) == 0:
					logger.warning('Could not find file for pattern: ' + pattern)
					continue

				#Download files which match the pattern
				s3_downloader(client, bucket, pattern_files, logger, force)

#--------------------------------------------------------------------------------------------------------#
def run_downloaddata(args):
	
	logger = TobiasLogger("DownloadData", args.verbosity)
	logger.begin()

	#Create config dict from commandline
	config = {"endpoint": args.endpoint, "username": args.username, "accesskey": args.key, "buckets": {args.bucket: args.patterns}}

	#Overwrite if args.yaml is set:
	if args.yaml is not None:
		config_from_yaml = read_config_yaml(args.yaml)
		for key in config:
			config[key] = config_from_yaml[key]
	
	logger.debug("Configuration dict is: {0}".format(config))

	#Download data using s3 client
	s3_client(config, logger, args.force)

	logger.end()
