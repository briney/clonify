#!/usr/local/bin/python
# filename: nr.py

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




import json
import os
import sqlite3
import subprocess as sp
import sys
import tempfile
import time




def make_nr_seq_db(args):
	db_path = os.path.join(args.temp, 'nr_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute('CREATE TABLE seqs (nr_seq_id text, all_seq_ids text)')
	return c


def make_lookup_db(args):
	print('\nMaking a JSON lookup database...')
	db_path = os.path.join(args.temp, 'lookup_db')
	conn = sqlite3.connect(db_path)
	c = conn.cursor()
	c.execute('DROP TABLE IF EXISTS seqs')
	c.execute('CREATE TABLE seqs (seq_id text, json blob)')
	return c


def fill_nr_seq_db(seq_db, seq_id_clusters):
	'''
	Fills a non-redundant sequence database.
	seq_db: a sqlite database connection object
	nr_seqs: a list of tuples, of the following format:
			 (representative_seq_id, [list of seq ids of identical sequences])
	'''
	print('Filling the non-redundant sequence database...')
	insert_data = [(s[0], ' '.join(s)) for s in seq_id_clusters]
	insert_cmd = 'INSERT INTO seqs VALUES (?,?)'
	seq_db.executemany(insert_cmd, insert_data)


def index_nr_db(nr_db):
	print('\nIndexing the non-redundant sequence database...')
	nr_db.execute('CREATE INDEX seq_index ON seqs (nr_seq_id)')


def fill_lookup_db(lookup_db, seqs):
	print('Filling the JSON lookup database...')
	pseqs = []
	for s in seqs:
		json_string = json.dumps(s)
		pseqs.append((s['seq_id'], json_string))
	insert_cmd = 'INSERT INTO seqs VALUES (?,?)'
	lookup_db.executemany(insert_cmd, pseqs)
	lookup_db.execute('CREATE INDEX seq_index ON seqs (seq_id)')


def chunker(l, size=900):
	return (l[pos:pos + size] for pos in range(0, len(l), size))


def json_lookup(seq_ids, lookup_db):
	print('Looking up non-redundant JSONs')
	results = []
	for chunk in chunker(seq_ids):
		result_chunk = lookup_db.execute('''SELECT seqs.json
								   FROM seqs
								   WHERE seqs.seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
		results.extend([r[0] for r in result_chunk.fetchall()])
	jsons = [json.loads(r) for r in results]
	return jsons


def expand_nr_seqs(seq_ids, nr_db):
	results = []
	for chunk in chunker(seq_ids):
		result_chunk = nr_db.execute('''SELECT seqs.all_seq_ids
								   FROM seqs
								   WHERE seqs.nr_seq_id IN ({})'''.format(','.join('?' * len(chunk))), chunk)
		flat_result_chunk = [r for sublist in result_chunk.fetchall() for r in sublist]
		results.extend(flat_result_chunk)
	seqs = [r.split() for r in results]
	return [s for sublist in seqs for s in sublist]


def remove_sqlite_db(name, args):
	db_path = os.path.join(args.temp, name)
	os.unlink(db_path)


def remove_nr_db(args):
	remove_sqlite_db('nr_db', args)


def remove_lookup_db(args):
	remove_sqlite_db('lookup_db', args)





# =========================================
#
#           CD-HIT CLUSTERING
#
# =========================================



def cdhit_clustering(seqs, args):
	infile = make_cdhit_input(seqs, args)
	outfile = os.path.join(args.temp, 'clust')
	logfile = open(os.path.join(args.temp, 'log'), 'a')
	do_cdhit(infile.name, outfile, logfile, args)
	clust_handle = open('{}.clstr'.format(outfile), 'r')
	cluster_ids = parse_clusters(clust_handle)
	if not args.debug:
		os.unlink(infile.name)
		os.unlink(os.path.join(args.temp, 'log'))
		os.unlink(outfile)
		os.unlink(outfile + '.clstr')
	return cluster_ids


def make_cdhit_input(seqs, args):
	infile = tempfile.NamedTemporaryFile(dir=args.temp, delete=False, mode='w')
	fastas = []
	for s in seqs:
		seq_type = 'vdj_nt' if args.non_redundant == 'nt' else 'vdj_aa'
		fastas.append('>{}\n{}'.format(s['seq_id'], s[seq_type]))
	infile.write('\n'.join(fastas))
	infile.close()
	return infile


def do_cdhit(fasta, clust, log, args):
	seq_type = 'nucleotide' if args.non_redundant == 'nt' else 'amino acid'
	print('\nClustering {} sequences with CD-HIT...'.format(seq_type), end='')
	sys.stdout.flush()
	start_time = time.time()
	cdhit_cmd = 'cd-hit -i {} -o {} -c {} -n 5 -d 0 -T 0 -M 35000'.format(
		fasta, clust, args.clustering_threshold)
	cluster = sp.Popen(cdhit_cmd, shell=True, stdout=log)
	cluster.communicate()
	print('Done.\nClustering took {} seconds'.format(round(time.time() - start_time, 2)))


def parse_clusters(cluster_handle):
	print('\nParsing CD-HIT cluster file...', end='')
	sys.stdout.flush()
	clusters = [c.split('\n') for c in cluster_handle.read().split('\n>')]
	print('Done.\n{} total clusters identified.\n'.format(len(clusters)))
	print('Retrieving cluster sequence IDs...', end='')
	sys.stdout.flush()
	all_cluster_ids = []
	start = time.time()
	lengths = []
	for cluster in clusters:
		lengths.append(len(cluster) - 1)
		cluster_ids = get_cluster_ids(cluster)
		all_cluster_ids.append(cluster_ids)
	print('Done.\n{} non-redundant clusters were identified)'.format(len(lengths)))
	print('The average cluster contains {} sequences; the largest contains {} sequences.'.format(round(1. * sum(lengths) / len(lengths), 2), max(lengths)))
	print('Cluster parsing and sequence retrieval took {} seconds\n'.format(round(time.time() - start, 2)))
	return all_cluster_ids


def get_cluster_ids(cluster):
	ids = []
	for c in cluster[1:]:
		if c:
			ids.append(c.split()[2][1:-3])
	return ids


def expand_clusters(nr_clusters, nr_db):
	index_nr_db(nr_db)
	print('Expanding non-redundant sequences...', end='')
	sys.stdout.flush()
	clusters = {}
	for cluster_id in list(nr_clusters.keys()):
		cluster_seq_ids = nr_clusters[cluster_id]
		expanded_ids = expand_nr_seqs(cluster_seq_ids, nr_db)
		clusters[cluster_id] = expanded_ids
	print('Done.')
	return clusters


def make_nr(seqs, nr_db, args):
	lookup_db = make_lookup_db(args)
	fill_lookup_db(lookup_db, seqs)
	seq_cluster_ids = cdhit_clustering(seqs, args)
	fill_nr_seq_db(nr_db, seq_cluster_ids)
	nr_ids = [s[0] for s in seq_cluster_ids]
	nr_jsons = json_lookup(nr_ids, lookup_db)
	remove_lookup_db(args)
	return nr_jsons
