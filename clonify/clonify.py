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

import os
import sys
import json
import time
import urllib
import sqlite3
import argparse
import subprocess as sp
import multiprocessing as mp
from collections import OrderedDict

from pymongo import MongoClient



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
args = parser.parse_args()



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


def get_db():
	if args.user and args.password:
		password = urllib.quote_plus(password)
		uri = 'mongodb://{}:{}@{}'.format(args.user, password, args.ip)
		conn = MongoClient(uri)
	else:
		conn = MongoClient(args.ip, 27017)
	return conn[args.db]


def get_collections():
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


def query(collection):
	c = db[collection]
	junc_query = 'junc_aa'
	vdj_query = 'vdj_nt' if args.non_redundant == 'nt' else 'vdj_aa'
	results = c.find({'chain': 'heavy', 'prod': 'yes'},
					 {'_id': 0, 'seq_id': 1, 'v_gene.full': 1, 'j_gene.full': 1, junc_query: 1, vdj_query: 1, 'var_muts_nt': 1})
	return [r for r in results]


def ensure_index(field, group):
	print('\nEnsuring indexes prior to updating:')
	for collection in group:
		print("Indexing '{}' on {}...".format(field, collection))
		coll = db[collection]
		coll.ensure_index(field)


def update_db(clusters, group):
	print_update_info()
	start = time.time()
	clust_count = len(clusters)
	p = mp.Pool(processes=250)
	async_results = []
	for clust_id in clusters:
		seq_ids = clusters[clust_id]
		async_results.append(p.apply_async(update, args=(clust_id, seq_ids, group)))
	monitor_update(async_results)
	p.close()
	p.join()
	run_time = time.time() - start
	print('Updating took {} seconds. ({} sequences per second)'.format(round(run_time, 2), round(len(clusters) / run_time, 1)))


def update(clust_id, seq_ids, group):
	clust_size = len(seq_ids)
	for collection in group:
		c = db[collection]
		c.update({'seq_id': {'$in': seq_ids}},
				 {'$set': {'clonify': {'id': clust_id, 'size': clust_size}}},
				 multi=True)





################################
#
#           SQLITE
#
################################



def build_sqlite_db(sequences):
	print('Building a SQLite database of sequences...', end='')
	sys.stdout.flush()
	seqs = [(s['seq_id'], s['junc_aa']) for s in sequences]
	db_path = os.path.join(args.temp, 'seq_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	create_cmd = get_seq_db_creation_cmd()
	insert_cmd = get_seq_db_insert_cmd()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute(create_cmd)
	c.executemany(insert_cmd, seqs)
	print('Done')
	print('Indexing the SQLite database...', end='')
	sys.stdout.flush()
	start = time.time()
	c.execute('CREATE INDEX seq_index ON seqs (seq_id)')
	print('Done')
	print('Indexing took {} seconds\n'.format(round(time.time() - start, 2)))
	return c


def get_seq_db_creation_cmd():
	return '''CREATE TABLE seqs (seq_id text, junction text)'''


def get_seq_db_insert_cmd():
	return 'INSERT INTO seqs VALUES (?,?)'


def remove_sqlite_db():
	db_path = os.path.join(args.temp, 'seq_db')
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






################################
#
#           CLONIFY
#
################################



def clonify(json_file):
	cluster_file = os.path.join(args.temp, '{}_clusters.txt'.format(args.db))
	clonify_cmd = './utils/cluster {} {}'.format(json_file, cluster_file)
	p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
	stdout, stderr = p.communicate()
	print(stdout.strip())
	if args.debug:
		print(stderr)
	clusters = []
	with open(cluster_file) as f:
		for line in f:
			clusters.append(line.strip())
	if not args.debug:
		os.unlink(json_file)
		os.unlink(cluster_file)
	return clusters


def build_clonify_input(seqs, json_file):
	'''
	Builds a JSON file suitable for input into Clonify.

	Input
	A pymongo query result object

	Output
	Path to a JSON-formatted file
	'''
	json_handle = open(json_file, 'a')
	jsons = []
	seq_ids = []
	for s in seqs:
		try:
			seq_id = s['seq_id']
			seq_ids.append(seq_id)
			v_full = s['v_gene']['full']
			v_fam = v_full.split('-')[0][4:]
			v_gene = '-'.join(v_full.split('*')[0].split('-')[1:])
			v_all = v_full.split('*')[1]
			j_full = s['j_gene']['full']
			j_fam = j_full.split('*')[0][4:]
			j_all = j_full.split('*')[1]
			junc_aa = s['junc_aa']
			pre_muts = s['var_muts_nt']['muts']
			muts = []
			for m in pre_muts:
				muts.append({'loc': str(m['loc']), 'mut': m['mut']})
			mut_count = s['var_muts_nt']['num']
			j = OrderedDict([
				('v_gene', {'all': v_all,
							'gene': v_gene,
							'full': v_full,
							'fam': v_fam}),
				('seq_id', seq_id),
				('j_gene', {'all': j_all,
							'full': j_full,
							'gene': j_fam}),
				('junc_aa', junc_aa),
				('var_muts_nt', {'muts': muts,
								 'num': mut_count})
			])
			jsons.append(j)
		except:
			continue
	json_handle.write(json.dumps(jsons, indent=2) + '\n')
	json_handle.close()
	return seq_ids


def parse_cluster_file(clust_ids, seq_ids):
	clusters = {}
	clust_ids = [c for c in clust_ids if c]
	for seq_id, clust_id in zip(seq_ids, clust_ids):
		lineage = clust_id
		if lineage in clusters:
			clusters[lineage].append(seq_id)
		else:
			clusters[lineage] = [seq_id, ]
	return {c: seqs for c, seqs in clusters.iteritems() if len(seqs) >= 2}


def write_output(clusters, sql, collection_group=None):
	collection = '_'.join(collection_group) if collection_group is not None else args.collection
	oname = '{}_{}_clusters.txt'.format(args.db, collection) if collection else '{}_clusters.txt'.format(args.db)
	ofile = os.path.join(args.output, oname)
	open(ofile, 'w').write('')
	ohandle = open(ofile, 'a')
	for c in clusters:
		clust_id = 'lineage_{}'.format(c)
		seqs = get_cluster_seqs(clusters[c], sql)
		ostring = '#{}\n{}\n\n'.format(clust_id, '\n'.join(seqs))
		ohandle.write(ostring)


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
	# 	if args.collection in get_collections():
	# 		print('\nFound the requested collection: {}'.format(args.collection))
	# else:
	# 	if args.collection_prefix:
	# 		print('\nSearching {} for collections that start with "{}"'.format(args.db, args.collection_prefix))
	# 	else:
	# 		print('\nSearching {} for collections'.format(args.db))
	# 	collections = get_collections()
	# 	print('\nFound {} collections:\n{}'.format(len(collections), ', '.join(collections)))
	print('\nThe following groups of collections will be processed:')
	print('\n'.join([', '.join(c) for c in coll_groups]))
	if args.non_redundant:
		print('\nIdentical sequences from each collection will be collapsed before lineage assignment.')
	if not args.output:
		print('\nWARNING: No output directory given so no output will be written.')
	if not args.update:
		print('\nWARNING: The MongoDB database will not be updated with clonality information.')
	print('\n\n')
	print('=' * 32)
	print('     PROCESSING COLLECTIONS')
	print('=' * 32)


def print_collection_info(collection):
	coll_len = len(collection)
	dash_len = len(collection) + 12
	print('')
	print('')
	print('-' * dash_len)
	print('     ' + collection)
	print('-' * dash_len)
	print('')
	print('Querying for sequences...')
	sys.stdout.flush()


def print_collection_group_info(group, num):
	gname_string = 'Collection Group #{}'.format(num)
	print('\n\n')
	print(gname_string)
	print('=' * len(gname_string))
	print('')
	print(', '.join(group))


def print_query_info(coll_seqs, nr=None):
	print('Found {} sequences'.format(len(coll_seqs)))
	if nr:
		print('{} were left after collapsing identical sequences.'.format(len(nr)))


def print_clonify_input_building(seqs):
	rstring = 'RUNNING CLONIFY'
	print('\n\n')
	print('=' * (len(rstring) + 16))
	print(rstring)
	print('=' * (len(rstring) + 16))
	print('\nBuilding Clonify input for {} total sequences.'.format(len(seqs)))


def print_clonify_start():
	print('\nRunning Clonify...')
	sys.stdout.flush()


def print_clonify_end(clusters):
	print('{} total clusters were identified.'.format(len(clusters)))


def print_output():
	print('\nWriting clonality information to file ({}).'.format(args.output))


def print_update_info():
	print('')
	print('Updating the MongoDB database with clonality info')


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
	print('\n{0} sequences were segregated into {1} clonal families.'.format(seq_count, clust_count))
	print('The average cluster size was %0.2f.' % (avg_clust_size))
	print('The largest cluster contains {} seqeunces.'.format(max_clust_size))
	print('%s sequences were assigned to clonal families (%0.2f%% of all sequences).' % (clustered_seqs, 100.0 * clustered_seqs / seq_count))
	print('\n')





def main():
	# print_start_info()
	# seq_ids = []
	# seq_dict = {}
	# json_file = os.path.join(args.temp, '{}.json'.format(args.db))


	collection_groups = get_collection_groups()
	print_start_info(collection_groups)

	for i, collection_group in enumerate(collection_groups):
		print_collection_group_info(collection_group, i)
		seqs = []
		json_file = os.path.join(args.temp, '{}_{}.json'.format(args.db, i))
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
		if len(seqs) == 0:
			continue
		print_clonify_input_building(seqs)
		seq_ids = build_clonify_input(seqs, json_file)
		print_clonify_start()
		cluster_ids = clonify(json_file)
		clusters = parse_cluster_file(cluster_ids, seq_ids)
		print_clonify_end(clusters)
		if args.non_redundant:
			clusters = nr.expand_clusters(clusters, nr_db)
			nr.remove_nr_db(args)
		if args.update:
			ensure_index('seq_id', collection_group)
			update_db(clusters, collection_group)
		if args.output:
			print_output()
			sql = build_sqlite_db(seqs)
			name = collection_group if len(collection_group) == 1 else str(i)
			write_output(clusters, sql, collection_group=name)
			remove_sqlite_db()
		print_finished(seq_ids, clusters)




	# if args.non_redundant:
	# 	from utils import nr
	# 	nr_db = nr.make_nr_seq_db(args)
	# 	for collection in get_collections():
	# 		print_collection_info(collection)
	# 		coll_seqs = query(collection)
	# 		nr_coll_seqs = nr.make_nr(coll_seqs, nr_db, args)
	# 		if args.pool:
	# 			if 'pool' in seq_dict:
	# 				seq_dict['pool'].extend(nr_coll_seqs)
	# 			else:
	# 				seq_dict['pool'] = nr_coll_seqs
	# 		else:
	# 			seq_dict[collection] = nr_coll_seqs
	# 		print_query_info(coll_seqs, nr=nr_coll_seqs)
	# else:
	# 	for collection in get_collections():
	# 		print_collection_info(collection)
	# 		coll_seqs = query(collection)
	# 		if args.pool:
	# 			if 'pool' in seq_dict:
	# 				seq_dict['pool'].extend(coll_seqs)
	# 			else:
	# 				seq_dict['pool'] = coll_seqs
	# 		else:
	# 			seq_dict[collection] = coll_seqs
	# 		print_query_info(coll_seqs)
	# for name in sorted(seq_dict.keys()):
	# 	seqs = seq_dict[name]
	# 	print_clonify_input_building(name, seqs)
	# 	seq_ids = build_clonify_input(seqs, json_file)
	# 	print_clonify_start()
	# 	cluster_ids = clonify(json_file)
	# 	clusters = parse_cluster_file(cluster_ids, seq_ids)
	# 	print_clonify_end(clusters)
	# 	if args.non_redundant:
	# 		clusters = nr.expand_clusters(clusters, nr_db)
	# 		nr.remove_nr_db(args)
	# 	if args.update:
	# 		ensure_index('seq_id')
	# 		update_db(clusters)
	# 	if args.output:
	# 		print_output()
	# 		sql = build_sqlite_db(seqs)
	# 		write_output(clusters, sql)
	# 	print_finished(seq_ids, clusters)



if __name__ == '__main__':
	db = get_db()
	main()
