from __future__ import print_function
import sys
import re
import datetime
import os, hashlib, base64

import logzero
from logzero import logger

DEFAULT_FORMAT = '%(color)s[%(levelname)1.1s %(asctime)s %(module)s:%(lineno)d]%(end_color)s %(message)s'
logfmt = logzero.LogFormatter(fmt=DEFAULT_FORMAT, datefmt="%Y-%m-%d %H:%M:%S")
logzero.setup_default_logger(formatter=logfmt)

class InputStream(object):
    '''This class handles opening either stdin or a gzipped or non-gzipped file'''

    def __init__(self, string):
        '''Create a new wrapper around a stream'''
        self.temp_path_to_delete = None
        if string in (None, '-', 'stdin') and self.valid(string):
            self.handle = sys.stdin
            return
        if string.startswith('gs:'):
            default_workspace = os.path.join(os.environ['HOME'], 'bulk-download')
            workspace = default_workspace
            import google.auth
            from google.cloud import storage
            #Note: this will only work on the cloud
            #If you have to run outside the cloud you could authenticate
            #with `gcloud auth application-default login` but this is
            #actually not a good idea to do as yourself.  Use a service
            #account if you have to do this.
            #See: https://cloud.google.com/sdk/gcloud/reference/auth/application-default/login
            credentials, project = google.auth.default()
            storage_client = storage.Client(credentials=credentials, project=project_id)
            if not os.path.exists(default_workspace):
                logger.info("Creating directory: {}".format(workspace))
                os.makedirs(default_workspace)
            (localpath, cloud_md5) = download_blob(string, google_storage_client, workspace)
            verify_download(localpath, cloud_md5)
            string = localpath
            self.temp_path_to_delete = localpath
        if string.endswith('.gz'):
            import gzip
            self.handle = gzip.open(string, 'rb')
        else:
            self.handle = open(string, 'r')

    def verify_download(localpath, cloud_md5):
        logger.info("Verifying download")
        hexdigest = md5(localpath)
        if (hexdigest != cloud_md5):
            msg = "Corrupted download: {} -- expected MD5: {} | observed MD5: {}"
            msg = msg.format(localpath, cloud_md5, hexdigest)
            logger.error(msg)
            sys.exit("[err] : {}".format(msg))

    def derive_local_path(cloud_file, workspace):
        basefile = os.path.basename(cloud_file)
        dstpath = os.path.join(workspace, basefile)
        return dstpath

    def download_blob(cloudpath, google_storage_client, workspace):
        bucket_name = os.path.dirname(cloudpath).split('/')[2]
        source_blob_name = '/'.join(cloudpath.split('/')[3:])
        dstpath = derive_local_path(cloudpath, workspace)
        bucket = google_storage_client.get_bucket(bucket_name)
        blob = bucket.get_blob(source_blob_name)
        if blob is None:
            msg = "Did not find blob: {} in bucket: {}"
            msg = msg.format(source_blob_name, bucket_name)
            logger.error(msg)
            sys.exit("[err] : {}".format(msg))
        logger.info("Start download blob to file system: {}".format(dstpath))
        blob.download_to_filename(dstpath)
        logger.info("Finished download blob to file system")
        cloud_md5 = base64.b64decode(blob.md5_hash).encode('hex')
        return (dstpath, cloud_md5)


    def readline(self):
        l = self.handle.readline()
        return l

    def cleanup(path):
        if os.path.exists(path):
            logger.info("Deleting File: {}".format(path))
            os.unlink(path)

    @staticmethod
    def valid(string):
        '''Check if we allow the opening of stdin'''
        if string is None and sys.stdin.isatty():
            raise IOError('no input specified but terminal is interactive')
        return True

    def __enter__(self):
        '''Support use of with by passing back the originating handle'''
        return self.handle

    def __exit__(self, *kwargs):
        '''Support use of with by closing on exit of the context'''
        self.handle.close()
        if self.temp_path_to_delete != None:
            cleanup(self.temp_path_to_delete)

    def __iter__(self):
        '''Support use in loops like a normal file object'''
        return self.handle.__iter__()

    def close(self):
        '''Close the underlying handle'''
        return self.handle.close()

def parse_bnd_alt_string(alt_string):
    '''
    Parse the BND alt string and return separators and region
    '''
    # NOTE The below is ugly but intended to match things like [2:222[ and capture the brackets
    result = re.findall(r'([][])(.+?)([][])', alt_string)
    assert result, "%s\n" % alt_string
    #sys.stderr.write("%s\n" % alt_string)
    sep1, region, sep2 = result[0]
    assert sep1 == sep2
    chrom2, breakpoint2 = region.rsplit(':', 1)
    breakpoint2 = breakpoint2
    return sep1, chrom2, breakpoint2
