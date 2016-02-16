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


from __future__ import print_function

import argparse
import celery
from collections import OrderedDict
import cPickle as pickle
import json
import multiprocessing as mp
import os
import sqlite3
import subprocess as sp
import sys
import tempfile
import time
import traceback
import urllib

from abtools import mongodb
from abtools.queue.celery import celery
from abtools.utils import progbar

from utils import cluster
from utils.cluster import Cluster
from utils.database import Database



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
    parser.add_argument('-t', '--temp', dest='temp', default=None,
                        help="The directory in which temp files will be stored. \
                        If the directory doesn't exist, it will be created. Required.")
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
    parser.add_argument('-x', '--dist', dest='distance_cutoff', default=0.26, type=float,
                        help="The cutoff adjusted edit distance (aED) for segregating sequences into clonal families. \
                        Defaults to 0.35.")
    parser.add_argument('-n', '--no_update', dest='update', action='store_false', default=True,
                        help="Use to skip updating the MongoDB database with clonality info.")
    parser.add_argument('-D', '--debug', dest='debug', action='store_true', default=False,
                        help="If set, will run in debug mode.")
    return parser.parse_args()


class Args(object):
    """docstring for Args"""
    def __init__(self, db=None, collection=None,
        collection_prefix=None, collection_prefix_split=None, collection_prefix_split_pos=0,
        split_num=1, pool=False, ip='localhost', port=27017, user=None, password=None,
        output='', temp=None, non_redundant=False, clustering_threshold=1.0,
        distance_cutoff=0.35, update=True, debug=False):
        super(Args, self).__init__()
        if any([db is None, temp is None]):
            print('ERROR: the following options are required:')
            print('--db, --temp')
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
        self.non_redundant = non_redundant
        self.clustering_threshold = clustering_threshold
        self.distance_cutoff = float(distance_cutoff)
        self.update = update
        self.debug = debug


################################
#
#            MONGO
#
################################


def get_collection_groups():
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


# def get_db():
#     if args.user and args.password:
#         password = urllib.quote_plus(password)
#         uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
#         conn = MongoClient(uri)
#     else:
#         conn = MongoClient(args.ip, 27017)
#     return conn[args.db]


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
    vdj_query = 'vdj_nt' if args.non_redundant == 'nt' else 'vdj_aa'
    results = c.find({'chain': 'heavy', 'prod': 'yes'},
                     {'_id': 0, 'seq_id': 1, 'v_gene.full': 1, 'j_gene.full': 1, 'junc_aa': 1, vdj_query: 1, 'var_muts_nt': 1})
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
    clust_count = len(clusters)
    p = mp.Pool(processes=250)
    async_results = []
    for c in clusters:
        async_results.append(p.apply_async(update, args=(c, group)))
    monitor_update(async_results)
    p.close()
    p.join()
    run_time = time.time() - start
    logger.info('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2), round(len(clusters) / run_time, 1)))

    # print_update_info()
    # start = time.time()
    # clust_count = len(clusters)
    # p = mp.Pool(processes=250)
    # async_results = []
    # for clust_id in clusters:
    #     seq_ids = clusters[clust_id]
    #     async_results.append(p.apply_async(update, args=(clust_id, seq_ids, group)))
    # monitor_update(async_results)
    # p.close()
    # p.join()
    # run_time = time.time() - start
    # logger.info('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2), round(len(clusters) / run_time, 1)))


def update(clust, group):
    for collection in group:
        c = db[collection]
        c.update({'seq_id': {'$in': clust.seq_ids}},
                 {'$set': {'clonify': {'id': clust.name, 'size': clust.size}}},
                 multi=True)

# def update(clust_id, seq_ids, group):
#     clust_size = len(seq_ids)
#     for collection in group:
#         c = db[collection]
#         c.update({'seq_id': {'$in': seq_ids}},
#                  {'$set': {'clonify': {'id': clust_id, 'size': clust_size}}},
#                  multi=True)


def get_sequences(collection_group, args):
    seqs = []
    if args.non_redundant:
        from utils import nr
        nr_db = nr.make_nr_seq_db(args)
        for collection in collection_group:
            print_collection_info(collection)
            coll_seqs = query(collection)
            nr_coll_seqs = nr.make_nr(coll_seqs, nr_db, args)
            seqs += nr_coll_seqs
            print_query_info(coll_seqs, nr=nr_coll_seqs)
    else:
        for collection in collection_group:
            print_collection_info(collection)
            coll_seqs = query(collection)
            seqs += coll_seqs
            print_query_info(coll_seqs)
    return seqs


################################
#
#           SQLITE
#
################################


def build_sqlite_db(sequences):
    logger.info('Building a SQLite database of sequences...')
    seqs = [(s['seq_id'], s['junc_aa']) for s in sequences]
    db_path = os.path.join(args.temp, 'seq_db')
    conn = sqlite3.connect(db_path)
    c = conn.cursor()
    create_cmd = get_seq_db_creation_cmd()
    insert_cmd = get_seq_db_insert_cmd()
    c.execute('DROP TABLE IF EXISTS seqs')
    c.execute(create_cmd)
    c.executemany(insert_cmd, seqs)
    build_sqlite_index(c, 'seqs', 'seq_id', 'seq_index')
    # logger.info('Indexing the SQLite database...')
    # start = time.time()
    # c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
    # logger.info('Indexing took {} seconds\n'.format(round(time.time() - start, 2)))
    return c


def build_json_db(seq_ids, json_files):
    logger.info('Building a SQLite database of JSON files...')
    jsondb = Database('json_db', args.temp_dir)
    data = [(j, pickle.dumps(s)) for j, s in zip(json_files, seq_ids)]
    jsondb.insert_many(data)
    jsondb.index()
    return jsondb


def build_mr_db(sequences):
    logger.info('Building a SQLite database of sequences for map-reduce lookup...')
    # hackety hackerson to pickle the JSON data into a SQLite db, but
    # since the db is just for rapid lookup of the JSON data by seq_id,
    # I'm OK with it.
    mrdb = Database('mr_db', args.temp_dir)
    seqs = [(s['seq_id'], pickle.dumps(s, protocol=0)) for s in sequences]
    mrdb.insert_many(seqs)
    mrdb.index()
    return mrdb

    # seqs = [(s['seq_id'], pickle.dumps(s, protocol=0)) for s in sequences]
    # db_path = os.path.join(args.temp, 'mr_db')
    # conn = sqlite3.connect(db_path)
    # c = conn.cursor()
    # create_cmd = get_seq_db_creation_cmd()
    # insert_cmd = get_seq_db_insert_cmd()
    # c.execute('DROP TABLE IF EXISTS seqs')
    # c.execute(create_cmd)
    # c.executemany(insert_cmd, seqs)
    # build_sqlite_index(c, 'seqs', 'seq_id', 'seq_index')
    # return c


def get_seq_db_creation_cmd():
    return '''CREATE TABLE seqs (seq_id text, junction text)'''


def get_seq_db_insert_cmd():
    return 'INSERT INTO seqs VALUES (?,?)'


def get_mr_db_creation_cmd():
    return '''CREATE TABLE seqs (seq_id text, json text)'''


def get_mr_db_insert_cmd():
    return 'INSERT INTO seqs VALUES (?,?)'


def build_sqlite_index(cursor, table, field, index_name):
    logger.info('Indexing the SQLite database...')
    start = time.time()
    cursor.execute('CREATE INDEX {} ON {} ({})'.format(index_name, table, field))
    logger.info('Indexing took {} seconds\n'.format(round(time.time() - start, 2)))


def remove_sqlite_db(name='seq_db'):
    db_path = os.path.join(args.temp, name)
    os.unlink(db_path)


def chunker(l, size=900):
    return (l[pos:pos + size] for pos in xrange(0, len(l), size))


def get_cluster_seqs(seq_ids, sql):
    seqs = []
    for chunk in chunker(seq_ids):
        seq_chunk = sql.execute('''SELECT seqs.seq_id, seqs.junction
                                   FROM seqs
                                   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
        seqs.extend(seq_chunk)
    return ['>{}\n{}'.format(s[0], s[1]) for s in seqs]


def get_mr_seqs(seq_ids, sql):
    seqs = []
    for chunk in chunker(seq_ids):
        seq_chunk = sql.execute('''SELECT seqs.seq_id, seqs.json
                                   FROM seqs
                                   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
        seqs.extend(seq_chunk)
    return [pickle.loads(str(s[1])) for s in seqs]


################################
#
#           CLONIFY
#
################################


def clonify(json_files, mr_db, json_db, args):
    if args.celery:
        cluster_files = run_map_jobs_via_celery(json_files, json_db, args.debug)
    else:
        cluster_files = run_map_jobs_via_multiprocessing(json_files, json_db, args.debug)
    return run_reduce_jobs_via_multiprocessing(cluster_files, mr_db, args)


def run_map_jobs_via_multiprocessing(json_files, json_db, debug):
    p = mp.Pool(maxtasksperchild=50)
    async_results = []
    for j in json_files:
        async_results.append(p.apply_async(clonify_map, (j, json_db, debug)))
    monitor_mp_jobs([a[1] for a in async_results])
    results = []
    for a in async_results:
        try:
            results.append(a[1].get())
        except:
            logger.debug('FILE-LEVEL EXCEPTION: {}'.format(a[0]))
            logging.debug(''.join(traceback.format_exc()))
    p.close()
    p.join()
    return results


def run_reduce_jobs_via_multiprocessing(cluster_files, mr_db, args):
    starting_file_count = len(cluster_files)
    while len(cluster_files) > 1:
        reduced_cluster_files = []
        cluster_sets = chunker(cluster_files, 2)
        for cs in cluster_sets:
            if len(cs) == 1:
                reduced_cluster_files += cs
            async_results.append(p.apply_async(clonify_reduce, (cs, mr_db)))
        reduced_cluster_files += [ar.get() for ar in async_results]
        cluster_files = reduced_cluster_files
    return cluster.get_clusters_from_file(cluster_files[0])


def monitor_mp_jobs(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        ready = [ar for ar in results if ar.ready()]
        finished = len(ready)
        progbar.progress_bar(finished, jobs)
    sys.stdout.write('\n')


def run_map_jobs_via_celery(json_files, json_db, debug):
    async_results = []
    for f in files:
        async_results.append((f, clonify_map.delay(j, json_db, debug)))
    succeeded, failed = monitor_celery_jobs([ar[1] for ar in async_results])
    failed_files = [ar[0] for ar in async_results if ar[1].failed()]
    # failed_files = [f for i, f in enumerate(files) if async_results[i].failed()]
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


@celery.task
def clonify_map(json_file, json_db, debug=False):
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

    seq_ids = json_db.find(json_file, unpickle=True)
    cluster_file = json_file + '_cluster'

    # need to name for the clonify C++ program, and should put it
    # in a location on my $PATH so that I can call it directly.
    clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
    p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    stdout, stderr = p.communicate()
    logger.info(stdout.strip())
    logger.debug(stderr)
    cluster_names = []
    with open(cluster_file) as f:
        for line in f:
            cluster_names.append(line.strip())
    clusters = ['\t'.join(z) for z in zip(seq_ids, cluster_names)]
    seq_cluster_file = json_file + '_seq-cluster'
    open(seq_cluster_file, 'w').write('\n'.join(clusters))
    if not debug:
        os.unlink(json_file)
        os.unlink(cluster_file)
    return seq_cluster_file
    # cluster_file = os.path.join(args.temp_dir, '{}_clusters.txt'.format(os.path.basename(json_file)))
    # # cluster_cmd_path = os.path.normpath('./utils/cluster')
    # clonify_cmd = './utils/cluster {} {}'.format(json_file, cluster_file)
    # p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
    # stdout, stderr = p.communicate()
    # logger.info(stdout.strip())
    # logger.debug(stderr)
    # clusters = []
    # with open(cluster_file) as f:
    #     for line in f:
    #         clusters.append(line.strip())
    # if not args.debug:
    #     os.unlink(json_file)
    #     os.unlink(cluster_file)
    # return clusters


@celery.task
def clonify_reduce(cluster_pair, mr_db, json_db):
    '''
    Given a pair of cluster file names, merge their clusters.

    Returns the name of a new cluster file containing the merged results
    and deletes the input cluster files.
    '''
    # get Cluster objects for each cluster file
    merged = cluster.get_clusters_from_file(cluster_pair[0], mr_db)
    comparison = cluster.get_clusters_from_file(cluster_pair[1], mr_db)
    # compare the clusters
    for cc in comparison:
        merge_happened = False
        for i, cm in enumerate(merged):
            json_file = cm.json(other=cc, as_file=True)
            cluster_file = clonify_map(json_file, json_db)
            _clusters = cluster.get_clusters_from_file(cluster_file)
            # If cm isn't in _clusters, some merging happened.
            # We need to remove cm from the merged list and replace
            # it with _clusters (which is likely to be just one cluster,
            # but may contain two new clusters)
            if cm not in _clusters:
                del merged[i]
                merged += _clusters
                merge_happened = True
                break
        # if we've checked all of the merged clusters and the comparison
        # cluster didn't merge with any of them, append to the merged list
        if not merge_happened:
            merged.append(cc)
    # set up a new, merged cluster file
    temp_dir = os.path.dirname(cluster_pair[0])
    new_cluster_file = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
    # remove the old cluster files
    for old_file in cluster_pair:
        os.unlink(old_file)
    # assemble/write the data for the merged cluster file
    cluster_assignments = []
    for i, m in enumerate(merged):
        assignments = zip(m.seq_ids, [str(i)] * len(seq_ids))
        cluster_assignments += ['\t'.join(a) for a in assignments]
    new_cluster_file.write('\n'.join(cluster_assignments))
    return new_cluster_file.name






# def build_clonify_input(seqs, json_file):
#     '''
#     Builds a JSON file suitable for input into Clonify.

#     Input
#     A pymongo query result object

#     Output
#     Path to a JSON-formatted file
#     '''
#     json_handle = open(json_file, 'a')
#     jsons = []
#     seq_ids = []
#     for s in seqs:
#         try:
#             seq_id = s['seq_id']
#             seq_ids.append(seq_id)
#             v_full = s['v_gene']['full']
#             v_fam = v_full.split('-')[0][4:]
#             v_gene = '-'.join(v_full.split('*')[0].split('-')[1:])
#             v_all = v_full.split('*')[1]
#             j_full = s['j_gene']['full']
#             j_fam = j_full.split('*')[0][4:]
#             j_all = j_full.split('*')[1]
#             junc_aa = s['junc_aa']
#             pre_muts = s['var_muts_nt']['muts']
#             muts = []
#             for m in pre_muts:
#                 muts.append({'loc': str(m['loc']), 'mut': m['mut']})
#             mut_count = s['var_muts_nt']['num']
#             j = OrderedDict([
#                 ('v_gene', {'all': v_all,
#                             'gene': v_gene,
#                             'full': v_full,
#                             'fam': v_fam}),
#                 ('seq_id', seq_id),
#                 ('j_gene', {'all': j_all,
#                             'full': j_full,
#                             'gene': j_fam}),
#                 ('junc_aa', junc_aa),
#                 ('var_muts_nt', {'muts': muts,
#                                  'num': mut_count})
#             ])
#             jsons.append(j)
#         except:
#             continue
#     json_handle.write(json.dumps(jsons, indent=2) + '\n')
#     json_handle.close()
#     return seq_ids


# def parse_cluster_file(clust_ids, seq_ids):
#     clusters = {}
#     clust_ids = [c for c in clust_ids if c]
#     for seq_id, clust_id in zip(seq_ids, clust_ids):
#         lineage = clust_id
#         if lineage in clusters:
#             clusters[lineage].append(seq_id)
#         else:
#             clusters[lineage] = [seq_id, ]
#     return {c: seqs for c, seqs in clusters.iteritems() if len(seqs) >= 2}


def write_output(clusters, mr_db, collection_group=None):
    collection = '|'.join(collection_group) if collection_group is not None else args.collection
    oname = '{}_{}_clusters.txt'.format(args.db, collection) if collection else '{}_clusters.txt'.format(args.db)
    ofile = os.path.join(args.output, oname)
    open(ofile, 'w').write('')
    ohandle = open(ofile, 'a')
    for c in clusters:
        seq_ids = c.seq_ids
        jsons = mrdb.find(seq_ids, unpickle=True)
        output = ['#' + c.name] + [j['junc_aa'] for j in jsons]
        ohandle.write('\n'.join(output))
        ohandle.write('\n\n')

    # collection = '_'.join(collection_group) if collection_group is not None else args.collection
    # oname = '{}_{}_clusters.txt'.format(args.db, collection) if collection else '{}_clusters.txt'.format(args.db)
    # ofile = os.path.join(args.output, oname)
    # open(ofile, 'w').write('')
    # ohandle = open(ofile, 'a')
    # for c in clusters:
    #     clust_id = 'lineage_{}'.format(c)
    #     seqs = get_cluster_seqs(clusters[c], sql)
    #     ostring = '#{}\n{}\n\n'.format(clust_id, '\n'.join(seqs))
    #     ohandle.write(ostring)


def monitor_update(results):
    finished = 0
    jobs = len(results)
    while finished < jobs:
        time.sleep(1)
        finished = len([r for r in results if r.ready()])
        update_progress(finished, jobs, sys.stdout)
    sys.stdout.write('\n\n')


def update_progress(finished, jobs, log, failed=None):
    pct = int(100. * finished / jobs)
    ticks = pct / 2
    spaces = 50 - ticks
    if failed:
        prog_bar = '\r({}/{}) |{}{}|  {}% ({}, {})'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct, finished - failed, failed)
    else:
        prog_bar = '\r({}/{}) |{}{}|  {}%'.format(finished, jobs, '|' * ticks, ' ' * spaces, pct)
    sys.stdout.write(prog_bar)
    sys.stdout.flush()





def print_start_info(coll_groups):
    # if args.collection:
    #     if args.collection in get_collections():
    #         print('\nFound the requested collection: {}'.format(args.collection))
    # else:
    #     if args.collection_prefix:
    #         print('\nSearching {} for collections that start with "{}"'.format(args.db, args.collection_prefix))
    #     else:
    #         print('\nSearching {} for collections'.format(args.db))
    #     collections = get_collections()
    #     print('\nFound {} collections:\n{}'.format(len(collections), ', '.join(collections)))
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


def print_clonify_input_building(seqs):
    rstring = 'RUNNING CLONIFY'
    logger.info('\n\n')
    logger.info('=' * (len(rstring) + 16))
    logger.info(rstring)
    logger.info('=' * (len(rstring) + 16))
    logger.info('\nBuilding Clonify input for {} total sequences.'.format(len(seqs)))


def print_clonify_start():
    logger.info('\nRunning Clonify...')
    sys.stdout.flush()


def print_clonify_end(clusters):
    logger.info('{} total clusters were identified.'.format(len(clusters)))


def print_output():
    logger.info('\nWriting clonality information to file ({}).'.format(args.output))


def print_update_info():
    logger.info('')
    logger.info('Updating the MongoDB database with clonality info')


def print_finished(seq_ids, clusters):
    seq_count = len(seq_ids)
    clust_count = len(clusters)
    clust_sizes = [len(clusters[c]) for c in clusters.keys()]
    clustered_seqs = sum(clust_sizes)
    if clust_count > 0:
        avg_clust_size = float(clustered_seqs) / clust_count
        max_clust_size = max(clust_sizes)
    else:
        avg_clust_size = 0
        max_clust_size = 0
    logger.info('\n{0} sequences were segregated into {1} clonal families.'.format(seq_count, clust_count))
    logger.info('The average cluster size was %0.2f.' % (avg_clust_size))
    logger.info('The largest cluster contains {} seqeunces.'.format(max_clust_size))
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


# def chunker(l, n):
#     """Yield successive n-sized chunks from l."""
#     for i in xrange(0, len(l), n):
#         yield l[i:i + n]


def main(args):
    global db
    db = mongodb.get_db(args.db, args.ip, args.port, args.user, args.password)
    collection_groups = get_collection_groups()
    print_start_info(collection_groups)

    for i, collection_group in enumerate(collection_groups):
        print_collection_group_info(collection_group, i)
        seqs = get_sequences(collection_group, args)
        mr_db = build_mr_db(seqs)
        # seqs = []
        # json_file = os.path.join(args.temp, '{}_{}.json'.format(args.db, i))
        # if args.non_redundant:
        #     from utils import nr
        #     nr_db = nr.make_nr_seq_db(args)
        #     for collection in collection_group:
        #         print_collection_info(collection)
        #         coll_seqs = query(collection)
        #         nr_coll_seqs = nr.make_nr(coll_seqs, nr_db, args)
        #         seqs += nr_coll_seqs
        #         print_query_info(coll_seqs, nr=nr_coll_seqs)
        # else:
        #     for collection in collection_group:
        #         print_collection_info(collection)
        #         coll_seqs = query(collection)
        #         seqs += coll_seqs
        #         print_query_info(coll_seqs)
        if len(seqs) == 0:
            continue

        # Get the number of cores -- if multiprocessing, just mp.cpu_count();
        # if Celery, need to get the number of worker processes
        if args.celery:
            app = celery.Celery()
            app.config_from_object('abtools.celeryconfig')
            stats = app.control.inspect().stats()
            cores = sum([v['pool']["max-concurrency"] for v in stats.values()])
        else:
            cores = mp.cpu_count()

        # Divide the input sequences into JSON subfiles
        json_files = []
        seq_ids = []
        for seq_chunk in chunker(seqs, cores):
            seq_ids.append([s['seq_id'] for s in seq_chunk])
            json_files.append(Cluster.pretty_json(seq_chunk,
                as_file=True,
                temp_dir=args.temp))
        json_db = build_json_db(seq_ids, json_files)

        # Run clonify jobs via multiprocessing or Celery
        # Should return a list of Cluster objects for each subfile
        clusters = clonify(json_files, mr_db, json_lookup, args)



        # print_clonify_input_building(seqs)
        # seq_ids = build_clonify_input(seqs, json_file)
        # print_clonify_start()
        # cluster_ids = clonify(json_file)
        # clusters = parse_cluster_file(cluster_ids, seq_ids)
        # print_clonify_end(clusters)
        # if args.non_redundant:
        #     clusters = nr.expand_clusters(clusters, nr_db)
        #     nr.remove_nr_db(args)
        # if args.update:
        #     ensure_index('seq_id', collection_group)
        #     update_db(clusters, collection_group)
        # if args.output is not None:
        #     print_output()
        #     sql = build_sqlite_db(seqs)
        #     name = collection_group if len(collection_group) == 1 else str(i)
        #     write_output(clusters, sql, collection_group=name)
        #     remove_sqlite_db()
        # print_finished(seq_ids, clusters)



def run_standalone(args):
    global logger
    logfile = get_logfile(args)
    log.setup_logging(logfile)
    logger = log.get_logger('clonify')
    main(args)


def run(**kwargs):
    args = Args(**kwargs)
    global logger
    log.get_logger('clonify')
    mai(args)


if __name__ == '__main__':
    args = parse_args()
    logfile = get_logfile(args)
    log.setup_logging(logfile)
    logger = log.get_logger('clonify')
    main(args)
