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
import celery
from collections import OrderedDict
import json
import math
import multiprocessing as mp
import os
import sqlite3
import subprocess as sp
import sys
import tempfile
import time
import traceback
import urllib.request, urllib.parse, urllib.error

from abutils.utils import log, mongodb, progbar
from abutils.utils.pipeline import make_dir
from abtools.queue.celery import celery

from .utils import cluster
from .utils.cluster import Cluster, Clusters
from .utils.database import Database

if sys.version_info[0] > 2:
    import pickle
else:
    import cPickle as pickle



def parse_args():
    parser = argparse.ArgumentParser("")
    parser.add_argument('-d', '--db', dest='db', required=True,
                        help="The MongoDB database to be queried. Required.")
    parser.add_argument('-c', '--collection', dest='collection', default=None,
                        help="The MongoDB collection to be queried. \
                        If ommitted, sequences from all collections in the database will be processed.")
    parser.add_argument('--collection-prefix', dest='collection_prefix', default=None,
                        help="If supplied, Clonify will process only collections beginning with <collection_prefix>.")
    parser.add_argument('--collection-prefix-split', default=None,
                        help="If supplied, will split all collection names at the <split-num> occurance of the supplied string \
                        and group all collections with identical prefixes. --pool=True is implied with this option.")
    parser.add_argument('--collection-prefix-split-pos', default=0, type=int,
                        help="If supplied, will group all collections that are identical for the first <collection-prefix-split-pos> \
                        characters. --pool=True is implied with this option.")
    parser.add_argument('--split-num', default=1, type=int,
                        help="With <collection-prefix-split>, collection names will be split at the <split-num> occurance of \
                        <collection-prefix-split>. Uses 1-based indexing. Default is 1.")
    parser.add_argument('--pool', dest='pool', default=False, action='store_true',
                        help="If set, all collections will be pooled and processed as a single group.")
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
    parser.add_argument('-l', '--log', dest='logfile',
                        help="Path to the log file. Required.")
    parser.add_argument('--non-redundant', default=False,
                        help="Collapses identical sequences prior to running Clonify. \
                        Stores redundant sequence info so that the complete redundant Mongo database will be \
                        updated with lineage info (non-redunant sequences are re-expanded prior to database update). \
                        Options are 'nt' or 'aa', which collapse identical nucleotide or amino acid sequences, respectively. \
                        Best for collections that contain redundant sequences \
                        and are large enough to cause clonify segfaults.")
    parser.add_argument('--clustering-threshold', default=1.0, type=float,
                        help="Clustering threshold to be used with the --non-redundant option. \
                        Default is 1.0.")
    # The following option ('-x') doesn't do anything at the moment.
    parser.add_argument('-x', '--dist', dest='distance_cutoff', default=0.35, type=float,
                        help="NOT YET IMPLEMENTED. The cutoff adjusted edit distance (aED) for segregating \
                        sequences into clonal families. Defaults to 0.35.")
    parser.add_argument('-C', '--celery', dest="celery", default=False, action='store_true',
                        help="NOT YET IMPLEMENTED. Use if performing computation on a Celery cluster. \
                        If set, input files will be split into many subfiles and passed \
                        to a Celery queue. If not set, input files will still be split, but \
                        will be distributed to local processors using multiprocessing.")
    parser.add_argument('-w', '--num-map-workers', dest='num_map_workers', type=int, default=1,
                        help='The number of map process that will be spawned. Default is 1. \
                        Set to 0 to use the max number of available cores, whether on a local \
                        machine using multiprocessing or on a Celery cluster.')
    parser.add_argument('-n', '--no_update', dest='update', action='store_false', default=True,
                        help="Use to skip updating the MongoDB database with clonality info.")
    parser.add_argument('--test-algo', action='store_true', default=False,
                        help='Tests whether the cluster program works. Useful for troubleshooting.')
    parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    return parser.parse_args()


class Args(object):
    """docstring for Args"""
    def __init__(self, db=None, collection=None,
        collection_prefix=None, collection_prefix_split=None, collection_prefix_split_pos=0,
        split_num=1, pool=False, ip='localhost', port=27017, user=None, password=None,
        output='', temp=None, log=None, non_redundant=False, clustering_threshold=1.0,
        distance_cutoff=0.35, celery=False, update=True, debug=False):
        super(Args, self).__init__()
        if any([db is None, temp is None, logfile is None]):
            print('ERROR: the following options are required:')
            print('--db, --temp, --log')
            sys.exit(1)
        self.db = db
        self.collection = collection
        self.collection_prefix = collection_prefix
        self.collection_prefix_split = collection_prefix_split
        self.collection_prefix_split_pos = int(collection_prefix_split_pos)
        self.split_num = int(split_num)
        self.pool = pool
        self.ip = ip
        self.port = int(port)
        self.user = user
        self.password = password
        self.output = output
        self.temp = temp
        self.logfile = log
        self.non_redundant = non_redundant
        self.clustering_threshold = clustering_threshold
        self.distance_cutoff = float(distance_cutoff)
        self.celery = celery
        self.update = update
        self.debug = debug


def validate_args(args):
    for d in [args.output, args.temp]:
        if d is not None:
            make_dir(d)


################################
#
#            MONGO
#
################################


def get_collection_groups(args):
    if args.collection:
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
    all_collections = db.collection_names(include_system_collections=False)
    if args.collection_prefix:
        prefix_collections = [c for c in all_collections if c.startswith(args.collection_prefix)]
        if args.pool:
            return [sorted(prefix_collections), ]
        return [[c, ] for c in sorted(prefix_collections)]
    if args.collection_prefix_split:
        args.pool = True
        s = args.collection_prefix_split
        prefix_collections = []
        prefixes = sorted(list(set(
            [s.join(c.split(s)[:args.split_num]) for c in all_collections])))
        for p in prefixes:
            prefix_collections.append(sorted([c for c in all_collections if c.startswith(p)]))
        return prefix_collections
    if args.collection_prefix_split_pos:
        args.pool = True
        pos = args.collection_prefix_split_pos
        prefix_collections = []
        prefixes = sorted(list(set([c[:pos] for c in all_collections])))
        for p in prefixes:
            prefix_collections.append(sorted([c for c in all_collections if c.startswith(p)]))
        return prefix_collections
    if args.pool:
        return [sorted(all_collections), ]
    return [[c, ] for c in sorted(all_collections)]


def get_collections(args):
    if args.collection:
        if os.path.isfile(args.collection):
            colls = []
            with open(args.collection) as f:
                for line in f:
                    colls.append(line.strip())
            return sorted(colls)
        else:
            return [args.collection, ]
    collections = db.collection_names(include_system_collections=False)
    if args.collection_prefix:
        collections = [c for c in collections if c.startswith(args.collection_prefix)]
    return sorted(collections)


def query(collection, args):
    c = db[collection]
    # vdj_query = 'vdj_nt' if args.non_redundant == 'nt' else 'vdj_aa'
    results = c.find({'chain': 'heavy', 'prod': 'yes', 'cdr3_len': {'$gte': 2}},
                     {'_id': 0, 'seq_id': 1, 'v_gene.full': 1, 'j_gene.full': 1, 'junc_aa': 1,
                      'vdj_nt': 1, 'vdj_aa': 1, 'var_muts_nt': 1})
    return [r for r in results]


def ensure_index(field, group):
    logger.info('\nEnsuring indexes prior to updating:')
    for collection in group:
        logger.info("Indexing '{}' on {}...".format(field, collection))
        coll = db[collection]
        coll.ensure_index(field)


def update_db(clusters, group):
    print_update_info()
    start = time.time()
    sizes = []
    p = mp.Pool(processes=250)
    async_results = []
    for c in clusters:
        sizes.append(c.size)
        async_results.append(p.apply_async(update, args=(c, group)))
    monitor_update(async_results)
    p.close()
    p.join()
    seq_count = sum(sizes)
    run_time = time.time() - start
    logger.info('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2), round(seq_count / run_time, 1)))
    return sizes


def update(clust, group):
    for collection in group:
        c = db[collection]
        c.update({'seq_id': {'$in': clust.seq_ids}},
                 {'$set': {'clonify': {'id': clust.name, 'size': clust.size}}},
                 multi=True)


def get_sequences(collection_group, args):
    seqs = []
    if args.non_redundant:
        from .utils import nr
        nr_db = nr.make_nr_seq_db(args)
        for collection in collection_group:
            print_collection_info(collection)
            coll_seqs = query(collection, args)
            nr_coll_seqs = nr.make_nr(coll_seqs, nr_db, args)
            seqs += nr_coll_seqs
            print_query_info(coll_seqs, nr=nr_coll_seqs)
    else:
        for collection in collection_group:
            print_collection_info(collection)
            coll_seqs = query(collection, args)
            seqs += coll_seqs
            print_query_info(coll_seqs)
    return seqs


################################
#
#           SQLITE
#
################################


def build_json_db(seq_ids, json_files, args):
    logger.info('')
    logger.info('Building a database of JSON files...')
    jsondb = Database('json_db', args.temp)
    jsondb.insert_many(list(zip(json_files, seq_ids)))
    # for j, s in zip(json_files, seq_ids):
    #     jsondb[j] = s
    jsondb.commit()
    logger.info('Indexing...')
    jsondb.index()
    jsondb.close()
    return jsondb


def build_mr_db(sequences, args):
    logger.info('')
    logger.info('Building a database of sequences...')
    # hackety hackerson to pickle the JSON data into a SQLite db, but
    # since the db is just for rapid lookup of the JSON data by seq_id,
    # I'm OK with it.
    mrdb = Database('mr_db', args.temp)
    seqs = [(s['seq_id'], s) for s in sequences]
    mrdb.insert_many(seqs)
    logger.info('Indexing...')
    mrdb.index()
    mrdb.close()
    return mrdb


def chunker(l, size=900):
    return (l[pos:pos + size] for pos in range(0, len(l), size))


################################
#
#           CLONIFY
#
################################


def clonify(json_files, mr_db, json_db, args):
    # map
    if args.celery:
        logger.info('')
        logger.info('Running Clonify jobs via Celery...')
        cluster_files = run_map_jobs_via_celery(json_files, json_db, args)
    elif any([args.debug, args.num_map_workers == 1]):
        logger.info('')
        logger.info('Running Clonify...')
        cluster_files = run_map_jobs_singlethreaded(json_files, json_db, args)
    else:
        logger.info('')
        logger.info('Running Clonify jobs via multiprocessing...')
        cluster_files = run_map_jobs_via_multiprocessing(json_files, json_db, args)
    logger.info('')
    # reduce
    if args.num_map_workers == 1:
        return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)
    else:
        reduced_clusters = clonify_reduce(cluster_files, mr_db)
        return reduced_clusters


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


def run_map_jobs_singlethreaded(json_files, json_db, args):
    results = []
    for i, j in enumerate(json_files):
        logger.debug(j)
        results.append(clonify_map(j, json_db, args))
        progbar.progress_bar(i + 1, len(json_files))
    return results


def run_reduce_jobs_singlethreaded(cluster_files, mr_db, json_db, args):
    starting_file_count = len(cluster_files)
    total_iterations = get_iteration_count(starting_file_count)
    update_progress(0, total_iterations, iteration=starting_file_count)
    iteration = 1
    while len(cluster_files) > 1:
        reduced_cluster_files = []
        logger.debug('CLUSTER FILES: {}'.format(', '.join(cluster_files)))
        cluster_sets = chunker(cluster_files, 2)
        for cs in cluster_sets:
            if len(cs) == 1:
                logger.debug('SINGLE CLUSTER FILE: {}'.format(cs[0]))
                reduced_cluster_files += cs
                continue
            logger.debug('CLUSTER SET:\n{}\n'.format('\n'.join(cs)))
            reduced_cluster_files.append(clonify_reduce(cs, mr_db, json_db, args))
        cluster_files = reduced_cluster_files
        update_progress(iteration, total_iterations, iteration=len(cluster_files))
        logger.debug('ITERATION {}: {} cluster files remaining'.format(iteration, len(cluster_files)))
        iteration += 1
    logger.info('')
    return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)


def run_map_jobs_via_multiprocessing(json_files, json_db, args):
    p = mp.Pool(maxtasksperchild=50)
    async_results = []
    for j in json_files:
        logger.debug(j)
        async_results.append([j, p.apply_async(clonify_map, (j, json_db, args))])
    monitor_mp_jobs([a[1] for a in async_results])
    results = []
    for a in async_results:
        try:
            results.append(a[1].get())
        except:
            logger.debug('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
            logger.debug(''.join(traceback.format_exc()))
    p.close()
    p.join()
    return results


def run_reduce_jobs_via_multiprocessing(cluster_files, mr_db, json_db, args):
    starting_file_count = len(cluster_files)
    total_iterations = get_iteration_count(starting_file_count)
    update_progress(0, total_iterations, iteration=starting_file_count)
    p = mp.Pool(maxtasksperchild=50)
    iteration = 1
    while len(cluster_files) > 1:
        async_results = []
        reduced_cluster_files = []
        cluster_sets = chunker(cluster_files, 2)
        for cs in cluster_sets:
            if len(cs) == 1:
                logger.debug('SINGLE CLUSTER FILE: {}'.format(cs[0]))
                reduced_cluster_files += cs
                continue
            logger.debug('CLUSTER SET:\n{}\n'.format('\n'.join(cs)))
            async_results.append(p.apply_async(clonify_reduce, (cs, mr_db, json_db, args)))
        reduced_cluster_files += [ar.get() for ar in async_results]
        cluster_files = reduced_cluster_files
        update_progress(iteration, total_iterations, iteration=len(cluster_files))
        logger.debug('ITERATION {}: {} cluster files remaining'.format(iteration, len(cluster_files)))
        iteration += 1
    logger.info('')
    return cluster.get_clusters_from_file(cluster_files[0], mr_db=mr_db)


def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        progbar.progress_bar(finished, jobs)
    sys.stdout.write('\n')


def run_map_jobs_via_celery(json_files, json_db, args):
    async_results = []
    for f in files:
        async_results.append((f, clonify_map.delay(j, json_db, args)))
    succeeded, failed = monitor_celery_jobs([ar[1] for ar in async_results])
    failed_files = [ar[0] for ar in async_results if ar[1].failed()]
    for ff in failed_files:
        logger.debug('FAILED FILE: {}'.format(f))
    return [s.get() for s in succeeded]


def monitor_celery_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        succeeded = [ar for ar in results if ar.successful()]
        failed = [ar for ar in results if ar.failed()]
        finished = len(succeeded) + len(failed)
        progbar.progress_bar(finished, jobs)
    sys.stdout.write('\n')
    return succeeded, failed


def get_iteration_count(starting_file_count):
    for i in range(starting_file_count):
        if 2**i >= starting_file_count:
            return i


@celery.task
def clonify_map(json_file, json_db, args):
    '''
    Runs clonify on the sequences contained in a JSON file.

    Returns the name of a new file, of the format:
        seq_id1    cluster_name1
        seq_id2    cluster_name2
        ...

    The new file will be named <json_file>_seq-cluster
    The intermediate file will be named <json_file>_cluster

    The JSON file and intermediate cluster file will be deleted
    unless debug == True.
    '''
    seq_ids = json_db[json_file]
    cluster_file = json_file + '_cluster'
    # need to name for the clonify C++ program, and should put it
    # in a location on my $PATH so that I can call it directly.
    clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
    p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.debug(stdout.strip())
    logger.debug(stderr)
    cluster_names = []
    with open(cluster_file) as f:
        for line in f:
            cluster_names.append(line.strip())
    clusters = ['\t'.join(z) for z in zip(seq_ids, cluster_names)]
    seq_cluster_file = json_file + '_seq-cluster'
    open(seq_cluster_file, 'w').write('\n'.join(clusters))
    # clean up
    if not args.debug:
        json_db.close()
        os.unlink(json_file)
        os.unlink(cluster_file)
    return seq_cluster_file


@celery.task
def clonify_reduce(cluster_files, mr_db):

    '''
    Given a pair of cluster file names, merge their clusters.

    Returns the name of a new cluster file containing the merged results
    and deletes the input cluster files.
    '''
    # build a Clusters object of clusters from the clonify_map output
    clusters = Clusters()
    for f in cluster_files:
        clusters.add(cluster.get_clusters_from_file(f, mr_db=mr_db))

    # configure the Clusters object
    logger.info('')
    logger.info('Getting cluster sequences...')
    clusters.get_sequences()
    logger.info('')
    logger.info('Getting cluster centroids...')
    clusters.get_centroids()
    logger.info('')
    logger.info('Building cluster database...')
    clusters.build_cluster_db()

    # reduce
    print('Reducing clusters...')
    reduced_clusters = clusters.reduce()
    return reduced_clusters


def write_output(clusters, mr_db, args, collection_group=None):
    collection = '|'.join(collection_group) if collection_group is not None else args.collection
    oname = '{}_{}_clusters.txt'.format(args.db, collection) if collection else '{}_clusters.txt'.format(args.db)
    ofile = os.path.join(args.output, oname)
    open(ofile, 'w').write('')
    ohandle = open(ofile, 'a')
    for c in clusters:
        sequences = mr_db.find(c.seq_ids)
        output = ['#' + c.name] + ['>{}\n{}'.format(s['seq_id'], s['junc_aa']) for s in sequences]
        ohandle.write('\n'.join(output))
        ohandle.write('\n\n')


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
    ticks = pct / 2
    spaces = 50 - ticks
    if iteration is not None:
        prog_bar = '\r({}/{}) |{}{}|  {}% (files left: {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, iteration)
    else:
        prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)
    sys.stdout.flush()


def print_start_info(coll_groups, args):
    logger.info('\nThe following groups of collections will be processed:')
    logger.info('\n'.join([', '.join(c) for c in coll_groups]))
    if args.non_redundant:
        logger.info('\nIdentical sequences from each collection will be collapsed before lineage assignment.')
    if not args.output:
        logger.info('\nWARNING: No output directory given so no output will be written.')
    if not args.update:
        logger.info('\nWARNING: The MongoDB database will not be updated with clonality information.')
    logger.info('\n\n')
    logger.info('=' * 32)
    logger.info('     PROCESSING COLLECTIONS')
    logger.info('=' * 32)


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
    logger.info('\nRunning Clonify...')
    sys.stdout.flush()


def print_clonify_end(clusters):
    logger.info('{} total clusters were identified.'.format(len(clusters)))


def print_output(args):
    logger.info('')
    logger.info('Writing clonality information to file ({}).'.format(args.output))


def print_update_info():
    logger.info('')
    logger.info('Updating the MongoDB database with clonality info')


def print_finished(clust_sizes):
    non_singletons = [c for c in clust_sizes if c > 1]
    seq_count = sum(clust_sizes)
    clust_count = len(non_singletons)
    clustered_seqs = sum(non_singletons)
    if clust_count > 0:
        avg_clust_size = float(clustered_seqs) / clust_count
        max_clust_size = max(non_singletons)
    else:
        avg_clust_size = 0
        max_clust_size = 0
    logger.info('\n{0} sequences were segregated into {1} clonal families.'.format(seq_count, clust_count))
    logger.info('The average cluster size was %0.2f.' % (avg_clust_size))
    logger.info('The largest cluster contains {} sequences.'.format(max_clust_size))
    logger.info('%s sequences were assigned to clonal families (%0.2f%% of all sequences).' % (clustered_seqs, 100.0 * clustered_seqs / seq_count))
    logger.info('\n')


def get_logfile(args):
    if args.logfile is None:
    	if args.output is None:
    		return os.path.normpath('./{}.log'.format(args.db))
    	else:
    		return os.path.normpath(os.path.join(args.output, './{}.log'.format(args.db)))
    else:
    	return os.path.normpath(args.logfile)


def main(args):
    global db
    db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
    collection_groups = get_collection_groups(args)
    print_start_info(collection_groups, args)

    for i, collection_group in enumerate(collection_groups):
        print_collection_group_info(collection_group, i)
        seqs = get_sequences(collection_group, args)
        seq_count = len(seqs)
        mr_db = build_mr_db(seqs, args)
        if len(seqs) == 0:
            continue

        # Get the number of cores -- if multiprocessing, just mp.cpu_count();
        # if Celery, need to get the number of worker processes
        if args.num_map_workers > 0:
            cores = args.num_map_workers
        elif args.celery:
            app = celery.Celery()
            app.config_from_object('abtools.celeryconfig')
            stats = app.control.inspect().stats()
            cores = sum([v['pool']["max-concurrency"] for v in list(stats.values())])
        else:
            cores = mp.cpu_count()
        logger.debug('CORES: {}'.format(cores))

        # Divide the input sequences into JSON subfiles
        chunksize = int(math.ceil(float(len(seqs)) / cores))
        logger.debug('CHUNKSIZE: {}'.format(chunksize))
        json_files = []
        seq_ids = []
        for seq_chunk in chunker(seqs, chunksize):
            seq_ids.append([s['seq_id'] for s in seq_chunk])
            json_files.append(Cluster.pretty_json(seq_chunk,
                as_file=True,
                temp_dir=args.temp))
        json_db = build_json_db(seq_ids, json_files, args)

        # Run clonify jobs via multiprocessing or Celery
        # Should return a list of Cluster objects for each subfile
        # if args.test_algo:
        #     test_clonify_algorithm(json_files, mr_db, json_db, args)
        # else:
        #     clusters = clonify(json_files, mr_db, json_db, args)
        print_clonify_input_building(seq_count)
        clusters = clonify(json_files, mr_db, json_db, args)


        if args.output:
            print_output(args)
            name = collection_group if len(collection_group) == 1 else str(i)
            write_output(clusters, mr_db, args, collection_group=name)
        if args.update:
            cluster_sizes = update_db(clusters, collection_group)
        else:
            cluster_sizes = [c.size for c in clusters]
        print_finished(cluster_sizes)


def run_standalone(args):
    validate_args(args)
    global logger
    logfile = get_logfile(args)
    log.setup_logging(logfile, debug=args.debug)
    logger = log.get_logger('clonify')
    main(args)


def run(**kwargs):
    validate_args(args)
    args = Args(**kwargs)
    global logger
    log.get_logger('clonify')
    mai(args)


if __name__ == '__main__':
    args = parse_args()
    validate_args(args)
    logfile = get_logfile(args)
    log.setup_logging(logfile, debug=args.debug)
    logger = log.get_logger('clonify')
    main(args)
