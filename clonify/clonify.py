#!/usr/local/bin/python
# filename: clonify.py

#
# Copyright (c) 2015 Bryan Briney
# License: The MIT license (http://opensource.org/licenses/MIT)
#
# Permission is hereby granted, free of charge, to any person obtaining a copy of this software
# and associated documentation files (the "Software"), to deal in the Software without restriction,
# including without limitation the rights to use, copy, modify, merge, publish, distribute,
# sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all copies or
# substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING
# BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
# NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
# DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
#




import argparse
# import celery
from collections import OrderedDict
from datetime import datetime
import itertools
import json
import math
import multiprocessing as mp
import os
import shutil
import sqlite3
import subprocess as sp
import sys
import tempfile
from threading import Thread
import time
import traceback
import urllib.request, urllib.parse, urllib.error

import numpy as np

from abutils.core.sequence import Sequence
from abutils.utils import log, mongodb, progbar
from abutils.utils.cluster import cluster
from abutils.utils.inputs import read_input
from abutils.utils.jobs import monitor_mp_jobs
from abutils.utils.pipeline import list_files, make_dir
# from abtools.queue.celery import celery

from .utils.lineage import Lineages, Lineage
# from .utils.nr import expand_nr_seqs, index_nr_db
# from .utils.cluster import Cluster, Clusters
from .utils.database import ClonifyDB

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
    import pickle
else:
    STR_TYPES = [str, unicode]
    import cPickle as pickle



def parse_args():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-j', '--json', dest='json', default=None,
                        help="Path to a JSON file or a directory containing JSON files.")
    parser.add_argument('-d', '--db', dest='db', default=None,
                        help="The MongoDB database to be queried.")
    parser.add_argument('-c', '--collection', dest='collection', default=None,
                        help="The MongoDB collection to be queried. \
                        If ommitted, sequences from all collections in the database will be processed.")
    parser.add_argument('--selection-prefix', dest='selection_prefix', default=None,
                        help="If supplied, Clonify will process only MongoDB collections or JSON files beginning with <selection_prefix>.")
    parser.add_argument('--selection-prefix-split', default=None,
                        help="If supplied, will split all collection or JSON file names at the <split-num> occurance of the supplied string \
                        and group all collections or JSON files with identical prefixes. --pool=True is implied with this option.")
    parser.add_argument('--selection-prefix-split-pos', default=0, type=int,
                        help="If supplied, will group all collections or JSON files that are identical for the first <selection-prefix-split-pos> \
                        characters. --pool=True is implied with this option.")
    parser.add_argument('--split-num', default=1, type=int,
                        help="With <selection-prefix-split>, collection or JSON file names will be split at the <split-num> occurance of \
                        <selection-prefix-split>. Uses 1-based indexing. Default is 1.")
    parser.add_argument('--pool', dest='pool', default=False, action='store_true',
                        help="If set, all collections or JSON files will be pooled and processed as a single group.")
    parser.add_argument('-i', '--ip', dest='ip', default='localhost',
                        help="The IP address of the MongoDB server. Defaults to 'localhost'.")
    parser.add_argument('-P', '--port', dest='port', default=27017, type=int,
                        help="The port used to connect to the MongoDB server. Defaults to '27017'.")
    parser.add_argument('-u', '--user', dest='user', default=None,
                        help="Username for the MongoDB server. Not used if not provided.")
    parser.add_argument('-p', '--password', dest='password', default=None,
                        help="Password for the MongoDB server. Not used if not provided.")
    parser.add_argument('-o', '--out', dest='output', default='',
                        help="Directory for the output files. Files will be named '<collection>_clones.txt'. \
                        Failing to provide an output directory will result in no output files being written.")
    parser.add_argument('-t', '--temp', dest='temp', default='/tmp',
                        help="The directory in which temp files will be stored. \
                        If the directory doesn't exist, it will be created. Default is '/tmp'.")
    parser.add_argument('-l', '--log', dest='logfile', default=None,
                        help="Path to the log file. Required.")
    # parser.add_argument('--non-redundant', default=False,
    #                     help="Collapses identical sequences prior to running Clonify. \
    #                     Stores redundant sequence info so that the complete redundant Mongo database will be \
    #                     updated with lineage info (non-redunant sequences are re-expanded prior to database update). \
    #                     Options are 'nt' or 'aa', which collapse identical nucleotide or amino acid sequences, respectively. \
    #                     Best for collections that contain redundant sequences \
    #                     and are large enough to cause clonify segfaults.")
    parser.add_argument('--no-preclustering', dest='preclustering', action='store_false', default=True,
                        help="If set, will not perform pre-clustering by V/J genes and CDR3 homology. \
                        For small datasets this may improve accuracy, but for large datasets this will result \
                        in much longer runtimes (and possibly failure).")
    parser.add_argument('--clustering-threshold', default=0.65, type=float,
                        help="Threshold to be used when clustering VJ groups of sequences. \
                        Default is 1.0.")
    parser.add_argument('--clustering-memory-allocation', default=800, type=int,
                        help='Amount of memory allocated to CD-HIT for clustering of VJ groups, in MB. Default is 800')
    # The following option ('-x') doesn't do anything at the moment.
    # parser.add_argument('-x', '--dist', dest='distance_cutoff', default=0.35, type=float,
    #                     help="NOT YET IMPLEMENTED. The cutoff adjusted edit distance (aED) for segregating \
    #                     sequences into clonal families. Defaults to 0.35.")
    parser.add_argument('-C', '--celery', dest="celery", default=False, action='store_true',
                        help="NOT YET IMPLEMENTED. Use if performing computation on a Celery cluster. \
                        If set, input files will be split into many subfiles and passed \
                        to a Celery queue. If not set, input files will still be split, but \
                        will be distributed to local processors using multiprocessing.")
    # parser.add_argument('-w', '--num-map-workers', dest='num_map_workers', type=int, default=1,
    #                     help='The number of map process that will be spawned. Default is 1. \
    #                     Set to 0 to use the max number of available cores, whether on a local \
    #                     machine using multiprocessing or on a Celery cluster.')
    parser.add_argument('-n', '--no_update', dest='update', action='store_false', default=True,
                        help="Use to skip updating the MongoDB database with clonality info.")
    # parser.add_argument('--test-algo', action='store_true', default=False,
    #                     help='Tests whether the cluster program works. Useful for troubleshooting.')
    parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    args = parser.parse_args()
    args_dict = vars(args)
    return Args(**args_dict)


class Args(object):
    """docstring for Args"""
    def __init__(self, json=None, sequences=None, db=None, collection=None,
        selection_prefix=None, selection_prefix_split=None, selection_prefix_split_pos=0,
        split_num=1, pool=False, ip='localhost', port=27017, user=None, password=None,
        output='', temp=None, logfile=None, non_redundant=False, clustering_threshold=0.65, preclustering=True,
        clustering_memory_allocation=800, distance_cutoff=0.35, celery=False, update=True, debug=False):
        
        super(Args, self).__init__()
        
        self.json = json
        self.sequences = sequences
        self.db = db
        self.collection = collection
        self.selection_prefix = selection_prefix
        self.selection_prefix_split = selection_prefix_split
        self.selection_prefix_split_pos = int(selection_prefix_split_pos)
        self.split_num = int(split_num)
        self.pool = pool
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.output = output
        self.temp = temp
        self.logfile = logfile
        # self.non_redundant = non_redundant
        self.clustering_threshold = clustering_threshold
        self.preclustering = preclustering
        self.clustering_memory_allocation = int(clustering_memory_allocation)
        # self.distance_cutoff = float(distance_cutoff)
        self.celery = celery
        self.update = update
        self.debug = debug


def validate_args(args):
    if all([args.json is None, args.sequences is None, args.db is None]):
        print('\n\n')
        print('ERROR: at least one of the following options is required:')
        print('--json, --db')
        print('or a list of Sequence objects must be provided through the Clonify API.')
        sys.exit(1)

    if any([args.temp is None, args.logfile is None]):
        print('\n\n')
        print('ERROR: the following options are required:')
        print('--temp, --log')
        sys.exit(1)

    for d in [args.output, args.temp]:
        if d is not None:
            make_dir(d)





################################
#
#            INPUTS
#
################################


def get_input_groups(args):
    if args.sequences is not None:
        return ['sequences', ]
    if args.db is not None:
        return get_collection_groups(args)
    else:
        return get_file_groups(args)


def get_sequences(group, args):
    if args.sequences is not None:
        return args.sequences
    if args.db is not None:
        seqs = read_input(args.db, 'mongodb', collection=group,
                          mongo_ip=args.ip, mongo_port=args.port,
                          mongo_user=args.user, mongo_password=args.password,
                          query=QUERY, projection=PROJECTION)
    elif args.json is not None:
        seqs = read_input(group, 'json')
    return seqs.as_generator





################################
#
#            JSON
#
################################


BASE_JSON = '''  {{
    "v_gene": {{
      "all": "{v_all}", 
      "gene": "{v_gene}", 
      "full": "{v_full}", 
      "fam": "{v_fam}"
    }}, 
    "seq_id": "{seq_id}", 
    "j_gene": {{
      "all": "{j_all}", 
      "full": "{j_full}", 
      "gene": "{j_gene}"
    }}, 
    "junc_aa": "{junc_aa}", 
    "var_muts_nt": {{
      "muts": [\n{mut_string}\n      ],
      "num": {mut_num}
    }}
  }}'''


BASE_MUT = '''        {{
          "loc": "{loc}", 
          "mut": "{mut}"
        }}'''


def get_file_groups(args):
    if os.path.isfile(args.json):
        return [[args.json], ]
    elif os.path.isdir(args.json):
        all_files = list_files(args.json)
        if args.selection_prefix is not None:
            files = [f for f in files if os.path.basename(f).startswith(args.selection_prefix)]
            if args.pool:
                return [sorted(files), ]
            else:
                return [[f, ] for f in sorted(files)]
        elif args.selection_prefix_split is not None:
            args.pool = True
            s = args.selection_prefix_split
            files = []
            prefixes = sorted(list(set([s.join(os.path.basename(f).split(s)[:args.split_num]) for f in all_files])))
            for p in prefixes:
                files.append(sorted([f for f in all_files if os.path.basename(f).startswith(p)]))
            return files
        elif args.selection_prefix_split_pos is not None:
            pos = args.selection_prefix_split_pos
            files = []
            prefixes = sorted(list(set([os.path.basename(f)[:pos] for f in all_files])))
            for p in prefixes:
                files.append(sorted([f for f in all_files if os.path.basename(f).startswith(p)]))
            return files
        elif args.pool:
            return [all_files, ]
        else:
            return [[f, ] for f in all_files]
    else:
        print('ERROR: The --json option accepts either a file or directory path.')
        print('You supplied: {}'.format(args.json))
        print('Please try again with a valid file/directory path.')
        sys.exit(1)


def pretty_json(sequences, as_file=False, temp_dir=None, raw=False, id_field='seq_id'):
    jsons = []
    for s in sequences:
        if len(s['var_muts_nt']['muts']) == 0:
            bases = ['A', 'C', 'G', 'T']
            s['var_muts_nt']['num'] = 1
            s['var_muts_nt']['muts'] = [{'loc': '{}.0'.format(random.randint(1, 300)),
                                        'mut': '{}>{}'.format(random.choice(bases), random.choice(bases))}]
        # this is to account for old/new versions of AbStar
        if 'loc' in s['var_muts_nt']['muts'][0]:
            mut_list = [BASE_MUT.format(loc=m['loc'], mut=m['mut']) for m in s['var_muts_nt']['muts']]
        else:
            mut_list = [BASE_MUT.format(loc=m['position'], mut='{}>{}'.format(m['was'], m['is'])) for m in s['var_muts_nt']['muts']]
        mut_string = ', \n'.join(mut_list)
        jsons.append(BASE_JSON.format(v_all=s['v_gene']['full'].split('*')[-1],
                                        v_gene='-'.join(s['v_gene']['full'].split('*')[0].split('-')[1:]),
                                        v_full=s['v_gene']['full'],
                                        v_fam=s['v_gene']['full'].split('-')[0].replace('IGHV', ''),
                                        seq_id=s[id_field],
                                        j_all=s['j_gene']['full'].split('*')[-1],
                                        j_full=s['j_gene']['full'],
                                        j_gene=s['j_gene']['full'].split('*')[0].replace('IGHJ', ''),
                                        junc_aa=s['junc_aa'],
                                        mut_string=mut_string,
                                        mut_num=s['var_muts_nt']['num']))
    if raw:
        return ', \n  '.join(jsons)
    json_string = '[\n  ' + ', \n  '.join(jsons) + '\n] \n'
    if as_file:
        temp_dir = temp_dir if temp_dir is not None else '/tmp'
        jfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False, mode='w')
        jfile.write(json_string)
        jfile.close()
        return jfile.name
    return json_string


def update_json(lineage_files, group, args):
    logger.info('')
    logger.info('Updating JSON files with clonality info')
    ldict = {}
    sdict = {}
    for lf in lineage_files:
        l = Lineage(lineage_file=lf)
        for seq_id in l.seq_ids:
            ldict[seq_id] = l.id
            sdict[seq_id] = l.size            
    # update the JSON files
    # TODO


    # for s in args.sequences:
    #     seq_id = s['seq_id']
    #     s['clonify'] = {'id': ldict[seq_id], 'size': sdict[seq_id]}





################################
#
#            MONGO
#
################################


QUERY = {'prod': 'yes', 'chain': 'heavy'}
PROJECTION = {'_id': 0, 'seq_id': 1, 'v_gene.gene': 1, 'j_gene.gene': 1, 'junc_aa': 1,
              'junc_nt': 1, 'vdj_nt': 1, 'vdj_aa': 1, 'var_muts_nt': 1}

def get_mongo_database(args):
    return mongodb.get_db(args.db, ip=args.ip, port=args.port,
                          user=args.user, password=args.password)


def close_mongo_database(db):
    db.client.close()


def get_collection_groups(args):
    if args.collection is not None:
        if os.path.isfile(args.collection):
            colls = []
            with open(args.collection) as f:
                for line in f:
                    colls.append(line.strip())
            if args.pool:
                return [sorted(colls), ]
            return [[c, ] for c in sorted(colls)]
        else:
            return [[args.collection, ], ]
    db = get_mongo_database(args)
    all_collections = mongodb.get_collections(db)
    if args.selection_prefix:
        prefix_collections = [c for c in all_collections if c.startswith(args.selection_prefix)]
        if args.pool:
            return [sorted(prefix_collections), ]
        return [[c, ] for c in sorted(prefix_collections)]
    if args.selection_prefix_split:
        args.pool = True
        s = args.selection_prefix_split
        prefix_collections = []
        prefixes = sorted(list(set(
            [s.join(c.split(s)[:args.split_num]) for c in all_collections])))
        for p in prefixes:
            prefix_collections.append(sorted([c for c in all_collections if c.startswith(p)]))
        return prefix_collections
    if args.selection_prefix_split_pos:
        args.pool = True
        pos = args.selection_prefix_split_pos
        prefix_collections = []
        prefixes = sorted(list(set([c[:pos] for c in all_collections])))
        for p in prefixes:
            prefix_collections.append(sorted([c for c in all_collections if c.startswith(p)]))
        return prefix_collections
    if args.pool:
        return [sorted(all_collections), ]
    return [[c, ] for c in sorted(all_collections)]


def ensure_index(db, field, group):
    logger.info('\nEnsuring indexes prior to updating:')
    for collection in group:
        logger.info("Indexing '{}' on {}...".format(field, collection))
        db[collection].ensure_index(field)


def update_mongodb(lineage_files, group, args):
    logger.info('')
    logger.info('Updating the MongoDB database with clonality info')
    db = get_mongo_database(args)
    ensure_index(db, 'seq_id', group)
    start = datetime.now()
    progbar.progress_bar(0, len(lineage_files), start_time=start)
    update_threads = 250
    for i in range(0, len(lineage_files), update_threads):
        tlist = []
        end = min([i + update_threads, len(lineage_files)])
        for lf in lineage_files[i:end]:
            l = Lineage(lineage_file=lf)
            t = Thread(target=mongoupdate, args=(db, l, group))
            t.start()
            tlist.append(t)
        for t in tlist:
            t.join()
        progbar.progress_bar(end, len(lineage_files))


def mongoupdate(db, lineage, group):
    for collection in group:
        c = db[collection]
        update_result = c.update_many({'seq_id': {'$in': lineage.seq_ids}},
                                      {'$set': {'clonify': {'id': lineage.id, 'size': lineage.size}}})
        try:
            debug = ['MONGODB UPDATE RESULTS', 'Group contents:']
            debug += group
            debug.append('matching records: {}'.format(update_result.matched_count))
            debug.append('modified records: {}'.format(update_result.modified_count))
            debug.append('lineage size: {}'.format(lineage.size))
            debug.append('')
            logger.debug('\n'.join(debug))
        except:
            logger.debug(traceback.format_exc())





################################
#
#          SEQUENCES
#
################################


def update_sequences(lineage_files, args):
    logger.info('')
    logger.info('Updating sequences with clonality info')
    ldict = {}
    sdict = {}
    for lf in lineage_files:
        l = Lineage(lineage_file=lf)
        for seq_id in l.seq_ids:
            ldict[seq_id] = l.id
            sdict[seq_id] = l.size            
    # update the sequence objects
    for s in args.sequences:
        seq_id = s['seq_id']
        s['clonify'] = {'id': ldict[seq_id], 'size': sdict[seq_id]}





################################
#
#           SQLITE
#
################################


def build_clonify_db(sequences, args):
    logger.info('')
    logger.info('Building a SQLite database of Sequence data...')
    clonifydb = ClonifyDB('clonify_db', args.temp)
    clonifydb.insert(sequences)
    clonifydb.commit()
    logger.info('Indexing...')
    logger.info('sequence id')
    clonifydb.index(fields='id')
    logger.info('vgene and jgene')
    clonifydb.index(fields=['vgene', 'jgene'])
    clonifydb.close()
    return clonifydb


def group_by_vj(clonify_db, args):
    # logger.info('')
    # logger.info('Grouping sequences by V/J gene use...')
    grouping_dir = os.path.join(args.temp, 'vj_groups')
    make_dir(grouping_dir)
    vs = clonify_db.find_distinct('vgene')
    js = clonify_db.find_distinct('jgene')
    start = datetime.now()
    total = len(vs) * len(js)
    count = 0
    for v, j in itertools.product(vs, js):
        progbar.progress_bar(count, total, start_time=start, extra_info='{}, {}    '.format(v, j))
        seqs = clonify_db.get_cdr3s_for_vj_group(v, j)
        if not seqs:
            count += 1
            continue
        gfile = os.path.join(grouping_dir, '{}_{}'.format(v, j))
        with open(gfile, 'w') as f:
            for s in seqs:
                f.write('>{}\n{}\n'.format(s[0], s[1]))
        count += 1
    progbar.progress_bar(total, total, start_time=start, completion_string='\n', complete=True, extra_info=' ' * 15)
    return grouping_dir





################################
#
#           CLONIFY
#
################################


def cluster_vj_groups(groups, clonify_db, args):
    '''
    Clusters groups of sequences using a CDR3 similarity threshold. 

    Args:

        groups (list): a list of file paths, each corresponding to a single VJ group 
                       (sequences using the same V/J genes)
        
        args (Args): runtime arguments


    Returns:

        list: a list of file paths, each corresponsing to a single cluster and containing all
              of the sequences (pickled Sequence objects) belonging to that cluster. Cluster files are
              named as follows: ``<V-gene>_<J-gene>_<uniqueID>
    '''

    ## TODO: make this work with Celery/multiprocessing, will need to limit clustering to a single thread,
    ## and may not produce massive speedups as large datasets clustered with a single thread may result in
    ## overall speed decreases
    ## Also, we remove the clustering_temp directory upon completion, which shouldn't be done if this
    ## function is going to be parallelized with Celery/multiprocessing

    cluster_dir = os.path.join(args.temp, 'vj_clusters')
    make_dir(cluster_dir)
    cluster_temp = os.path.join(args.temp, 'clustering_temp')
    make_dir(cluster_temp)
    start = datetime.now()
    for i, group in enumerate(groups):
        v, j = os.path.basename(group).split('_')
        progbar.progress_bar(i, len(groups), start_time=start, extra_info='{}, {}    '.format(v, j))
        clusters = cluster(group, threshold=args.clustering_threshold, return_just_seq_ids=True,
                           make_db=False, max_memory=args.clustering_memory_allocation, quiet=True)
        for num, id_list in enumerate(clusters):
            cluster_file = os.path.join(cluster_dir, '{}_{}_{}'.format(v, j, num))
            seqs = clonify_db.get_sequences_by_id(id_list)
            with open(cluster_file, 'wb') as f:
                pickle.dump(seqs, f, protocol=2)
    progbar.progress_bar(len(groups), len(groups), start_time=start, extra_info='{}, {}    '.format(v, j))
    if not args.debug:
        shutil.rmtree(cluster_temp)
    return cluster_dir


def clonify(seq_files, lineage_dir, args):
    '''
    runs clonify on a set of sequence files (each separately, as they represent clustered VJ groups).

    Runs either single-threaded, via multiprocessing or via Celery.
    '''
    if args.celery:
        logger.info('')
        logger.info('Running Clonify jobs via Celery...')
        sizes = run_clonify_via_celery(seq_files, lineage_dir, args)
    elif any([args.debug, ]):
        logger.info('')
        logger.info('Running Clonify...')
        sizes = run_clonify_singlethreaded(seq_files, lineage_dir, args)
    else:
        logger.info('')
        logger.info('Running Clonify jobs via multiprocessing...')
        sizes = run_clonify_via_multiprocessing(seq_files, lineage_dir, args)
    return lineage_dir, sizes


def run_clonify_singlethreaded(seq_files, lineage_dir, args):
    start = datetime.now()
    sizes = []
    fcount = len(seq_files)
    for i, seq_file in enumerate(seq_files):
        progbar.progress_bar(i, fcount, start_time=start)
        sizes += run_clonify(seq_file, lineage_dir, args)
    progbar.progress_bar(fcount, fcount, start_time=start,
                         complete=True, completion_string='\n')
    return sizes


def run_clonify_via_multiprocessing(seq_files, lineage_dir, args):
    start = datetime.now()
    p = mp.Pool(maxtasksperchild=50)
    async_results = []
    for f in seq_files:
        async_results.append([f, p.apply_async(run_clonify, (f, lineage_dir, args))])
    monitor_mp_jobs([a[1] for a in async_results])
    sizes = []
    for a in async_results:
        try:
            sizes += a[1].get()
        except:
            logger.debug('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
            logger.debug(''.join(traceback.format_exc()))
    p.close()
    p.join()
    return sizes


def run_clonify_via_celery(seq_files, lineage_dir, args):
    pass


def run_clonify(seq_file, lineage_dir, args):
    '''
    Runs clonify on the sequences contained in a JSON file.

    Returns a Clonify Lineages object.

    The JSON file and intermediate cluster file will be deleted
    unless ``args.debug == True``.
    '''
    with open(seq_file, 'rb') as f:
        seqs = pickle.load(f)
    seq_ids = [s['seq_id'] for s in seqs]
    if len(seq_ids) == 1:
        return []
    clonify_dir = os.path.join(args.temp, 'clonify')
    make_dir(clonify_dir)
    json_file = pretty_json(seqs, as_file=True, temp_dir=clonify_dir)
    cluster_file = json_file + '_cluster'
    # need to name for the clonify C++ program, and should put it
    # in a location on my $PATH so that I can call it directly.
    clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
    p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.debug(stdout.strip())
    logger.debug(stderr)
    lineages = Lineages(seq_ids, cluster_file)
    sizes = []
    # lineage_files = []
    for lineage in lineages:
        sizes.append(lineage.size)
        lineage.write(lineage_dir)
    # clean up
    if not args.debug:
        os.unlink(json_file)
        os.unlink(cluster_file)
    return sizes


def write_clonify_input(sequences, args):
    seq_dir = os.path.join(args.temp, 'sequences')
    make_dir(seq_dir)
    seq_file = os.path.join(seq_dir, 'sequences')
    with open(seq_file, 'wb') as f:
        pickle.dump(sequences, f)
    return seq_dir
     


def update_clonify_info(lineage_files, group, args):
    if args.db is not None:
        update_mongodb(lineage_files, group, args)
    elif args.json is not None:
        update_json(lineage_files, group, args)
    else:
        update_sequences(lineage_files, args)


def write_clonify_output(lineage_files, clonify_db, group_id, args):
    make_dir(args.output)
    ofile = os.path.join(args.output, 'group_{}'.format(group_id))
    with open(ofile, 'w') as f:
        for lf in lineage_files:
            l = Lineage(lineage_file=lf)
            seqs = clonify_db.get_sequences_by_id(l.seq_ids)
            fstring = '\n'.join(['>{}\n{}'.format(s['seq_id'], s['junc_aa']) for s in seqs])
            f.write('#{}\n'.format(l.id) + fstring + '\n\n')




# OLD VERSION 
# -------------

# def clonify(json_files, mr_db, json_db, args):
#     # map
#     if args.celery:
#         logger.info('')
#         logger.info('Running Clonify jobs via Celery...')
#         cluster_files = run_map_jobs_via_celery(json_files, json_db, args)
#     elif any([args.debug, args.num_map_workers == 1]):
#         logger.info('')
#         logger.info('Running Clonify...')
#         cluster_files = run_map_jobs_singlethreaded(json_files, json_db, args)
#     else:
#         logger.info('')
#         logger.info('Running Clonify jobs via multiprocessing...')
#         cluster_files = run_map_jobs_via_multiprocessing(json_files, json_db, args)
#     logger.info('')
#     # reduce
#     if args.num_map_workers == 1:
#         return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)
#     else:
#         reduced_clusters = clonify_reduce(cluster_files, mr_db)
#         return reduced_clusters


# def test_clonify_algorithm(json_files, mr_db, json_db, args):
#     p = mp.Pool(maxtasksperchild=50)
#     async_results = []
#     for j in json_files:
#         logger.debug(j)
#         async_results.append([j, p.apply_async(test_algo, (j, ))])
#     monitor_mp_jobs([a[1] for a in async_results])
#     results = []
#     for a in async_results:
#         try:
#             results.append(a[1].get())
#         except:
#             logger.info('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
#             logger.info(''.join(traceback.format_exc()))
#     p.close()
#     p.join()


# def test_algo(json_file):
#     cluster_file = json_file + '_cluster'
#     clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
#     p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
#     stdout, stderr = p.communicate()
#     try:
#         of = open(cluster_file, 'r')
#         of.close()
#     except IOError:
#         logger.info('\nFAILED FILE: {}'.format(json_file))
#         print(stderr.strip())
#         logger.debug(traceback.format_exc())


# def run_map_jobs_singlethreaded(json_files, json_db, args):
#     results = []
#     for i, j in enumerate(json_files):
#         logger.debug(j)
#         results.append(clonify_map(j, json_db, args))
#         progbar.progress_bar(i + 1, len(json_files))
#     return results


# def run_reduce_jobs_singlethreaded(cluster_files, mr_db, json_db, args):
#     starting_file_count = len(cluster_files)
#     total_iterations = get_iteration_count(starting_file_count)
#     update_progress(0, total_iterations, iteration=starting_file_count)
#     iteration = 1
#     while len(cluster_files) > 1:
#         reduced_cluster_files = []
#         logger.debug('CLUSTER FILES: {}'.format(', '.join(cluster_files)))
#         cluster_sets = chunker(cluster_files, 2)
#         for cs in cluster_sets:
#             if len(cs) == 1:
#                 logger.debug('SINGLE CLUSTER FILE: {}'.format(cs[0]))
#                 reduced_cluster_files += cs
#                 continue
#             logger.debug('CLUSTER SET:\n{}\n'.format('\n'.join(cs)))
#             reduced_cluster_files.append(clonify_reduce(cs, mr_db, json_db, args))
#         cluster_files = reduced_cluster_files
#         update_progress(iteration, total_iterations, iteration=len(cluster_files))
#         logger.debug('ITERATION {}: {} cluster files remaining'.format(iteration, len(cluster_files)))
#         iteration += 1
#     logger.info('')
#     return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)


# def run_map_jobs_via_multiprocessing(json_files, json_db, args):
#     p = mp.Pool(maxtasksperchild=50)
#     async_results = []
#     for j in json_files:
#         logger.debug(j)
#         async_results.append([j, p.apply_async(clonify_map, (j, json_db, args))])
#     monitor_mp_jobs([a[1] for a in async_results])
#     results = []
#     for a in async_results:
#         try:
#             results.append(a[1].get())
#         except:
#             logger.debug('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
#             logger.debug(''.join(traceback.format_exc()))
#     p.close()
#     p.join()
#     return results


# def run_reduce_jobs_via_multiprocessing(cluster_files, mr_db, json_db, args):
#     starting_file_count = len(cluster_files)
#     total_iterations = get_iteration_count(starting_file_count)
#     update_progress(0, total_iterations, iteration=starting_file_count)
#     p = mp.Pool(maxtasksperchild=50)
#     iteration = 1
#     while len(cluster_files) > 1:
#         async_results = []
#         reduced_cluster_files = []
#         cluster_sets = chunker(cluster_files, 2)
#         for cs in cluster_sets:
#             if len(cs) == 1:
#                 logger.debug('SINGLE CLUSTER FILE: {}'.format(cs[0]))
#                 reduced_cluster_files += cs
#                 continue
#             logger.debug('CLUSTER SET:\n{}\n'.format('\n'.join(cs)))
#             async_results.append(p.apply_async(clonify_reduce, (cs, mr_db, json_db, args)))
#         reduced_cluster_files += [ar.get() for ar in async_results]
#         cluster_files = reduced_cluster_files
#         update_progress(iteration, total_iterations, iteration=len(cluster_files))
#         logger.debug('ITERATION {}: {} cluster files remaining'.format(iteration, len(cluster_files)))
#         iteration += 1
#     logger.info('')
#     return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)


# def monitor_mp_jobs(results):
#     finished = 0
#     jobs = len(results)
#     while finished < jobs:
#         time.sleep(1)
#         ready = [ar for ar in results if ar.ready()]
#         finished = len(ready)
#         progbar.progress_bar(finished, jobs)
#     sys.stdout.write('\n')


# def run_map_jobs_via_celery(json_files, json_db, args):
#     async_results = []
#     for f in files:
#         async_results.append((f, clonify_map.delay(j, json_db, args)))
#     succeeded, failed = monitor_celery_jobs([ar[1] for ar in async_results])
#     failed_files = [ar[0] for ar in async_results if ar[1].failed()]
#     for ff in failed_files:
#         logger.debug('FAILED FILE: {}'.format(f))
#     return [s.get() for s in succeeded]


# def monitor_celery_jobs(results):
#     finished = 0
#     jobs = len(results)
#     while finished < jobs:
#         time.sleep(1)
#         succeeded = [ar for ar in results if ar.successful()]
#         failed = [ar for ar in results if ar.failed()]
#         finished = len(succeeded) + len(failed)
#         progbar.progress_bar(finished, jobs)
#     sys.stdout.write('\n')
#     return succeeded, failed


# def get_iteration_count(starting_file_count):
#     for i in range(starting_file_count):
#         if 2**i >= starting_file_count:
#             return i


# # @celery.task
# def clonify_map(json_file, json_db, args):
#     '''
#     Runs clonify on the sequences contained in a JSON file.

#     Returns the name of a new file, of the format:
#         seq_id1    cluster_name1
#         seq_id2    cluster_name2
#         ...

#     The new file will be named <json_file>_seq-cluster
#     The intermediate file will be named <json_file>_cluster

#     The JSON file and intermediate cluster file will be deleted
#     unless debug == True.
#     '''
#     seq_ids = json_db[json_file]
#     cluster_file = json_file + '_cluster'
#     # need to name for the clonify C++ program, and should put it
#     # in a location on my $PATH so that I can call it directly.
#     clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
#     p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
#     stdout, stderr = p.communicate()
#     logger.debug(stdout.strip())
#     logger.debug(stderr)
#     cluster_names = []
#     with open(cluster_file) as f:
#         for line in f:
#             cluster_names.append(line.strip())
#     clusters = ['\t'.join(z) for z in zip(seq_ids, cluster_names)]
#     seq_cluster_file = json_file + '_seq-cluster'
#     open(seq_cluster_file, 'w').write('\n'.join(clusters))
#     # clean up
#     if not args.debug:
#         json_db.close()
#         os.unlink(json_file)
#         os.unlink(cluster_file)
#     return seq_cluster_file


# # @celery.task
# def clonify_reduce(cluster_files, mr_db):

#     '''
#     Given a pair of cluster file names, merge their clusters.

#     Returns the name of a new cluster file containing the merged results
#     and deletes the input cluster files.
#     '''
#     # build a Clusters object of clusters from the clonify_map output
#     clusters = Clusters()
#     for f in cluster_files:
#         clusters.add(cluster.get_clusters_from_file(f, mr_db=mr_db))

#     # configure the Clusters object
#     logger.info('')
#     logger.info('Getting cluster sequences...')
#     clusters.get_sequences()
#     logger.info('')
#     logger.info('Getting cluster centroids...')
#     clusters.get_centroids()
#     logger.info('')
#     logger.info('Building cluster database...')
#     clusters.build_cluster_db()

#     # reduce
#     print('Reducing clusters...')
#     reduced_clusters = clusters.reduce()
#     return reduced_clusters


# def write_output(clusters, mr_db, args, collection_group=None):
#     collection = '|'.join(collection_group) if collection_group is not None else args.collection
#     oname = '{}_{}_clusters.txt'.format(args.db, collection) if collection else '{}_clusters.txt'.format(args.db)
#     ofile = os.path.join(args.output, oname)
#     open(ofile, 'w').write('')
#     ohandle = open(ofile, 'a')
#     for c in clusters:
#         sequences = mr_db.find(c.seq_ids)
#         output = ['#' + c.name] + ['>{}\n{}'.format(s['seq_id'], s['junc_aa']) for s in sequences]
#         ohandle.write('\n'.join(output))
#         ohandle.write('\n\n')


def monitor_update(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        finished = len([r for r in results if r.ready()])
        update_progress(finished, jobs)
    sys.stdout.write('\n\n')


def update_progress(finished, jobs, iteration=None):
    pct = int(100. * finished / jobs)
    ticks = int(pct / 2)
    spaces = 50 - ticks
    if iteration is not None:
        prog_bar = '\r({}/{}) |{}{}|  {}% (files left: {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, iteration)
    else:
        prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)
    sys.stdout.flush()


def print_start_info(groups, args):
    gtype = 'sequences'
    if args.db is not None:
        gtype = 'MongoDB collections'
    if args.json is not None:
        gtype = 'json files'
    logger.info('')
    logger.info('The following groups of {} will be processed:'.format(gtype))
    for group in groups:
        logger.info(', '.join(group))
    # if args.non_redundant:
    #     logger.info('\nIdentical sequences from each collection will be collapsed before lineage assignment.')
    if not args.output:
        logger.info('')
        logger.info('WARNING: No output directory given so no output will be written.')
    if not args.update:
        logger.info('')
        logger.info('WARNING: The MongoDB database will not be updated with clonality information.')
    logger.info('')
    logger.info('')
    header = '    PROCESSING {}    '.format(gtype.upper())
    logger.info('=' * len(header))
    logger.info(header)
    logger.info('=' * len(header))
    logger.info('')


def print_collection_info(collection):
    coll_len = len(collection)
    dash_len = len(collection) + 12
    logger.info('')
    logger.info('')
    logger.info('-' * dash_len)
    logger.info('     ' + collection)
    logger.info('-' * dash_len)
    logger.info('')
    logger.info('Querying for sequences...')
    sys.stdout.flush()


def print_collection_group_info(group, num):
    gname_string = 'Collection Group #{}'.format(num)
    logger.info('\n\n')
    logger.info(gname_string)
    logger.info('=' * len(gname_string))
    logger.info('')
    logger.info(', '.join(group))


def print_query_info(coll_seqs, nr=None):
    logger.info('Found {} sequences'.format(len(coll_seqs)))
    if nr:
        logger.info('{} were left after collapsing identical sequences.'.format(len(nr)))


def print_clonify_input_building(seq_count):
    rstring = 'RUNNING CLONIFY'
    logger.info('\n\n')
    logger.info('=' * (len(rstring) + 16))
    logger.info(rstring)
    logger.info('=' * (len(rstring) + 16))
    logger.info('\nBuilding Clonify input for {} total sequences.'.format(seq_count))


def print_clonify_start():
    # logger.info('\nRunning Clonify...')
    sys.stdout.flush()


def print_clonify_end(clusters):
    logger.info('{} total clusters were identified.'.format(len(clusters)))


def print_output(args):
    logger.info('')
    logger.info('Writing clonality information to file ({}).'.format(args.output))


# def print_update_info():
#     logger.info('')
#     logger.info('Updating the MongoDB database with clonality info')



# def print_finished(clust_sizes):
#     non_singletons = [c for c in clust_sizes if c > 1]
#     seq_count = sum(clust_sizes)
#     clust_count = len(non_singletons)
#     clustered_seqs = sum(non_singletons)
#     if clust_count > 0:
#         avg_clust_size = float(clustered_seqs) / clust_count
#         max_clust_size = max(non_singletons)
#     else:
#         avg_clust_size = 0
#         max_clust_size = 0
#     logger.info('\n{0} sequences were segregated into {1} clonal families.'.format(seq_count, clust_count))
#     logger.info('The average cluster size was %0.2f.' % (avg_clust_size))
#     logger.info('The largest cluster contains {} sequences.'.format(max_clust_size))
#     logger.info('%s sequences were assigned to clonal families (%0.2f%% of all sequences).' % (clustered_seqs, 100.0 * clustered_seqs / seq_count))
#     logger.info('\n')


def print_group_info(group, num, num_groups, args):
    '''
    Prints information about the currently in-process sequence/file/collection group.

    Args:

        grpup (list): list of sequences, JSON file paths, or MongoDB collectinos.

        num (int): the number of the current group.

        num_groups (int): total number of groups to be processed.

    Returns:

        None
    '''
    gtype = 'SEQUENCE'
    if args.db is not None:
        gtype = 'COLLECTION'
    if args.json is not None:
        gtype = 'JSON FILE'
    gstring = '  {} GROUP #{}  '.format(gtype, num)
    logger.info('=' * len(gstring))
    logger.info(gstring)
    logger.info('=' * len(gstring))
    logger.info('')
    logger.info('group {} of {}'.format())
    logger.info('')
    logger.info('\n'.join(group))
    logger.info('')


def print_clonify_results(seq_count, lineage_sizes):




    print(lineage_sizes)
    
    
    
    
    
    gt1_sizes = [l for l in lineage_sizes if l > 1] 
    lineage_count = len(gt1_sizes)
    mean = np.mean(gt1_sizes)
    assigned_seq_count = sum(gt1_sizes)
    assigned_seq_pct = 100.0 * assigned_seq_count / seq_count
    max_size = max(gt1_sizes)
    logger.info('{} sequences were segregated into {} clonal families.'.format(seq_count, lineage_count))
    logger.info('The mean cluster size was {:0.2f}.'.format(mean))
    logger.info('The largest cluster contains {} sequences.'.format(max_size))
    logger.info('{} sequences were assigned to clonal families ({:0.2f}% of all sequences).'.format(assigned_seq_count,
                                                                                                    assigned_seq_pct))


def get_logfile(args):
    if args.logfile is None:
    	if args.output is None:
    		return os.path.abspath('./{}.log'.format(args.db))
    	else:
    		return os.path.abspath(os.path.join(args.output, './{}.log'.format(args.db)))
    else:
    	return os.path.abspath(args.logfile)


def main(args):

    groups = get_input_groups(args)
    num_groups = len(groups)
    print_start_info(groups, args)

    lineage_dir = os.path.join(args.temp, 'lineages')
    make_dir(lineage_dir)

    for i, group in enumerate(groups, 1):

        if num_groups > 1:
            print_group_info(group, i, num_groups, args)

        logger.info('----------------')
        logger.info('  INITIALIZING  ')
        logger.info('----------------')
        logger.info('Getting sequences...')
        sequences = get_sequences(group, args)
        logger.info('Building a SQLite sequence database...')
        clonify_db = build_clonify_db(sequences, args)
        seq_count = clonify_db.count
        logger.info('Retrieved {} sequences'.format(seq_count))
        logger.info('')

        logger.info('------------------')
        logger.info('  PRE-CLUSTERING  ')
        logger.info('------------------')
        if args.preclustering:
            logger.info('Grouping sequences by V/J gene...')
            vj_group_dir = group_by_vj(clonify_db, args)
            vj_group_files = list_files(vj_group_dir)
            logger.info('Clustering sequences within each VJ group...')
            cluster_dir = cluster_vj_groups(vj_group_files, clonify_db, args)
            cluster_files = list_files(cluster_dir)
        else:
            logger.info('Clonify is being run without pre-clustering.')
            logger.info('Writing Clonify input...')
            cluster_dir = write_clonify_input(sequences, args)
            cluster_files = list_files(cluster_dir)
        logger.info('')

        logger.info('-----------')
        logger.info('  CLONIFY  ')
        logger.info('-----------')
        # logger.info('Running Clonify...')
        lineage_dir, lineage_sizes = clonify(cluster_files, lineage_dir, args)
        lineage_files = list_files(lineage_dir)
        print_clonify_results(seq_count, lineage_sizes)
        logger.info('')

        logger.info('----------')
        logger.info('  UPDATE  ')
        logger.info('----------')
        update_clonify_info(lineage_files, group, args)
        logger.info('')
        
        if args.update:
            logger.info('----------')
            logger.info('  OUTPUT  ')
            logger.info('----------')
            write_clonify_output(lineage_files, clonify_db, i, args)
            logger.info('')
        
        
        logger.info('---------------')
        logger.info('  CLEANING UP  ')
        logger.info('---------------')
        if not args.debug:
            logger.info('Removing VJ group files...')
            shutil.rmtree(vj_group_dir)
            logger.info('Removing VJ cluster files...')
            shutil.rmtree(cluster_dir)
            logger.info('Removing lineage files...')
            shutil.rmtree(lineage_dir)
            logger.info('Removing ClonifyDB...')
            clonify_db.destroy()
        else:
            logger.info('Debug mode is ON, so temporary files are not removed.')
    
        logger.info('')
        logger.info('')
        logger.info('')










    # #==============  OLD  ===============#
    
    # global db
    # db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
    # collection_groups = get_collection_groups(args)
    # print_start_info(collection_groups, args)

    # for i, collection_group in enumerate(collection_groups):
    #     print_collection_group_info(collection_group, i)
    #     seqs, nr_db = get_sequences(collection_group, args)
    #     seq_count = len(seqs)
    #     mr_db = build_mr_db(seqs, args)
    #     if len(seqs) == 0:
    #         continue

    #     # Get the number of cores -- if multiprocessing, just mp.cpu_count();
    #     # if Celery, need to get the number of worker processes
    #     if args.num_map_workers > 0:
    #         cores = args.num_map_workers
    #     elif args.celery:
    #         app = celery.Celery()
    #         app.config_from_object('abtools.celeryconfig')
    #         stats = app.control.inspect().stats()
    #         cores = sum([v['pool']["max-concurrency"] for v in list(stats.values())])
    #     else:
    #         cores = mp.cpu_count()
    #     logger.debug('CORES: {}'.format(cores))

    #     # Divide the input sequences into JSON subfiles
    #     chunksize = int(math.ceil(float(len(seqs)) / cores))
    #     logger.debug('CHUNKSIZE: {}'.format(chunksize))
    #     json_files = []
    #     seq_ids = []
    #     for seq_chunk in chunker(seqs, chunksize):
    #         seq_ids.append([s['seq_id'] for s in seq_chunk])
    #         json_files.append(Cluster.pretty_json(seq_chunk,
    #             as_file=True,
    #             temp_dir=args.temp))
    #     json_db = build_json_db(seq_ids, json_files, args)

    #     # Run clonify jobs via multiprocessing or Celery
    #     # Should return a list of Cluster objects for each subfile
    #     # if args.test_algo:
    #     #     test_clonify_algorithm(json_files, mr_db, json_db, args)
    #     # else:
    #     #     clusters = clonify(json_files, mr_db, json_db, args)
    #     print_clonify_input_building(seq_count)
    #     clusters = clonify(json_files, mr_db, json_db, args)

    #     if args.non_redundant:
    #         print('\nIndexing non-redundant sequence database...')
    #         index_nr_db(nr_db)
    #         print('Expanding clusters...')
    #         progbar.progress_bar(0, len(clusters))
    #         for i, c in enumerate(clusters, 1):
    #             expanded_seq_ids = expand_nr_seqs(c.seq_ids, nr_db)
    #             c.seq_ids = expanded_seq_ids
    #             c.size = len(expanded_seq_ids)
    #             progbar.progress_bar(i, len(clusters))
    #         print('\n')

    #     if args.output:
    #         print_output(args)
    #         name = collection_group if len(collection_group) == 1 else str(i)
    #         write_output(clusters, mr_db, args, collection_group=name)
    #     if args.update:
    #         ensure_index('seq_id', collection_group)
    #         cluster_sizes = update_db(clusters, collection_group)
    #     else:
    #         cluster_sizes = [c.size for c in clusters]
    #     print_finished(cluster_sizes)


def run_standalone(args):
    validate_args(args)
    global logger
    logfile = get_logfile(args)
    log.setup_logging(logfile, debug=args.debug)
    logger = log.get_logger('clonify')
    main(args)


def run(input=None, **kwargs):
    if input is not None:
        # check for a JSON file/directory
        if any([os.path.isfile(input), os.path.isdir(input)]):
            kwargs['json'] = input
        # check for an iterable of Sequence objects
        elif type(input) in [list, tuple] and type(input[0]) == Sequence:
            kwargs['sequences'] = input
        # only option left is a MongoDB directory
        else:
            kwargs['db'] = input
    args = Args(**kwargs)
    global logger
    log.get_logger('clonify')
    main(args)


if __name__ == '__main__':
    args = parse_args()
    validate_args(args)
    logfile = get_logfile(args)
    log.setup_logging(logfile, debug=args.debug)
    logger = log.get_logger('clonify')
    main(args)
