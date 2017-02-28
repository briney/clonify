#!/usr/local/bin/python
# filename: cluster.py

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


import collections
import itertools
import math
import os
import random
import string
import subprocess as sp
import tempfile

from database import Database

from abtools.cluster import cluster
from abtools.utils import progbar



class Clusters(object):
    """Manages a set of Clusters"""
    def __init__(self, clusters=None):
        super(Clusters, self).__init__()
        self._clusters = clusters
        self._db = None
        self._names = None
        self._centroids = None
        self._seq_ids = None


    @property
    def clusters(self):
        if self._clusters is None:
            return []
        return self._clusters

    @clusters.setter
    def clusters(self, clusters):
        self._clusters = clusters


    @property
    def centroids(self):
        if self._centroids is None:
            return []
        return self._centroids

    @centroids.setter
    def centroids(self, centroids):
        self._centroids = centroids


    @property
    def db(self):
        '''
        db is an in-memory SQLite key/value database with cluster names as
        keys and lists of sequence IDs as values.
        '''
        if self._db is None:
            self._db = Database(name='clusters', in_memory=True)
            self._db.index()
            # if self.clusters is not None:
            #     self.add_clusters_to_db(self.clusters)
        return self._db


    @property
    def names(self):
        return [c.name for c in self.clusters]


    def add(self, clusters):
        self.add_clusters_to_db(clusters)
        self.clusters += clusters


    def reduce(self):
        print('Writing reduce input file...')
        json_file = Cluster.pretty_json(self.centroids, as_file=True, id_field='name')
        seq_ids = [s['name'] for s in self.centroids]
        cluster_ids = self.clonify(json_file)
        print('Reducing clusters')
        reduced_clusters = {}
        for n, i in zip(cluster_ids, seq_ids):
            if n in reduced_clusters:
                reduced_clusters[n].add(i)
            else:
                reduced_clusters[n] = ReducedCluster(i)
        clusters = reduced_clusters.values()
        print('Getting reduced cluster sequence IDs...')
        progbar.progress_bar(0, len(clusters))
        for i, c in enumerate(clusters):
            ids = self.db.find(c.sub_clusters)
            ids = [_i for _i in itertools.chain(*ids)]
            # ids = []
            # for sc in c.sub_clusters:
            #     ids.extend(self.db.find_one(sc))
            # c.seq_ids = [_i for subl in ids for _i in subl]
            c.seq_ids = ids
            progbar.progress_bar(i + 1, len(clusters))
        return clusters



    def get_centroids(self):
        centroids = []
        total = len(self.clusters)
        progbar.progress_bar(0, total)
        for i, c in enumerate(self.clusters):
            centroids.extend(c.centroids)
            progbar.progress_bar(i + 1, total)
        self.centroids = centroids
        print('')


    def get_sequences(self):
        total = len(self.clusters)
        progbar.progress_bar(0, total)
        for i, c in enumerate(self.clusters):
            s = c.sequences
            progbar.progress_bar(i + 1, total)
        print('')


    def add_clusters_to_db(self, clusters):
        self.db.insert_many((c.name, c.seq_ids) for c in clusters)


    def build_cluster_db(self):
        if self._db is None:
            self._db = Database(name='clusters', in_memory=True)
            if self.clusters is not None:
                self.add_clusters_to_db(self.clusters)

    # def merge(self, other):
    #     cluster_names = self.names + other.names
    #     json_file = self._make_combined_clonify_input(other)
    #     cluster_ids = self._clonify(json_file)
    #     assignments = {}
    #     for n, i in zip(cluster_names, cluster_ids):
    #         assignments[n] = assignments[n].append(i) if n in assignments else [i, ]
    #     merged = [v for v in assignments.values() if len(v) > 1]
    #     _merged = False
    #     for m in merged:
    #         if cluster.name not in merged:
    #             continue
    #         _merged = True
    #         clusters = [cluster] + db.find(m)
    #         self._reconcile_merged_clusters(clusters)
    #     if not _merged:
    #         self.add_clusters_to_db([cluster])


    # def _reconcile_merged_clusters(self, clusters):
    #     sequences = []
    #     for c in clusters:
    #         sequences.extend(c.sequences)
    #     new_cluster = Cluster(sequences=sequences)
    #     self.db.delete([c.name for c in clusters])
    #     self.db.insert_one([new_cluster.name, new_cluster])
    #     self.clusters = list(self.db.find_all())



    # def _make_combined_clonify_input(self, other):
    #     sequences = self.centroids + other.centroids
    #     return Cluster.pretty_json(sequences, as_file=True)


    @staticmethod
    def clonify(json_file):
        cluster_file = json_file + '_cluster'
        clonify_cmd = 'cluster {} {}'.format(json_file, cluster_file)
        p = sp.Popen(clonify_cmd, stdout=sp.PIPE, stderr=sp.PIPE, shell=True)
        stdout, stderr = p.communicate()
        cluster_names = []
        with open(cluster_file) as f:
            for line in f:
                cluster_names.append(line.strip())
        # clean up
        os.unlink(json_file)
        os.unlink(cluster_file)
        return cluster_names




class Cluster(object):
    """docstring for Cluster"""
    def __init__(self, seq_ids=None, sequences=None, mrdb=None):
        super(Cluster, self).__init__()
        if seq_ids is not None:
            self.seq_ids = seq_ids
            self._sequences = None
        if sequences is not None:
            self._sequences = sequences
            self.seq_ids = [s['seq_id'] for s in sequences]
        self.db = mrdb
        self.size = len(self.seq_ids)
        self.json_file = None
        self._name = None
        self._v_gene = None
        self._cdr3_len = None
        self._centroids = None

    def __len__(self):
        return len(self.sequences)

    def __iter__(self):
        return iter(self.sequences)

    def __contains__(self, item):
        return item in self.seq_ids

    def __eq__(self, other):
        if hasattr(other, 'seq_ids'):
            return sorted(self.seq_ids) == sorted(other.seq_ids)
        return False


    @property
    def sequences(self):
        if self._sequences is None and self.db is not None:
            self._sequences = self.db.find(self.seq_ids)
        return self._sequences

    @sequences.setter
    def sequences(self, sequences):
        self._sequences = sequences

    @property
    def name(self):
        if self._name is None:
            self._name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        return self._name

    @property
    def v_gene(self):
        if self._v_gene is None:
            vgenes = [s['v_gene']['full'].split('*')[0] for s in self.sequences]
            vcounts = collections.Counter(vgenes)
            self._v_gene = sorted(vcounts.keys(), key=lambda x: vcounts[x], reverse=True)[0]
        return self._v_gene

    @property
    def cdr3_len(self):
        if self._cdr3_len is None:
            lengths = [len(s['junc_aa']) - 2 for s in self.sequences]
            lcounts = collections.Counter(lengths)
            self._cdr3_len = sorted(lcounts.keys(), key=lambda x: lcounts[x], reverse=True)[0]
        return self._cdr3_len


    @property
    def centroids(self):
        if self._centroids is None:
            # self._centroids = []
            # if self.size == 1:
            #     centroid = self.sequences[0]
            # else:
            #     c = cluster([(s['seq_id'], s['vdj_nt']) for s in self.sequences],
            #                 threshold=0.7,
            #                 quiet=True)
            #     c = sorted(c, key=lambda x: x.size, reverse=True)[0]
            #     cent = c.centroid
            #     centroid = self.db.find_one(cent.id)
            # centroid['name'] = '{}_{}'.format(self.name, '0')
            # self._centroids.append(centroid)

            num = int(math.ceil(self.size / 100.))
            for i in range(num):
                seq = random.choice(self.sequences)
                seq['name'] = '{}_{}'.format(self.name, i)
                self._centroids.append(seq)
        return self._centroids


    # def get_sequences(self):
    #     if self._sequences is None and self.db is not None:
    #         self._sequences = self.db.find(self.seq_ids)


    def should_compare(self, other):
        if all([self.v_gene == other.v_gene, abs(self.cdr3_len - other.cdr3_len) <= 5]):
            return True
        return False


    def json(self, other=None, as_file=False, temp_dir=None):
        other = other if other is not None else []
        sequences = self.sequences + other.sequences
        if as_file:
            self.json_file = self.pretty_json(sequences, as_file, temp_dir)
            return self.json_file
        else:
            return self.pretty_json(sequences, as_file, temp_dir)


    @staticmethod
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
            jfile = tempfile.NamedTemporaryFile(dir=temp_dir, delete=False)
            jfile.write(json_string)
            jfile.close()
            return jfile.name
        return json_string


    def delete_json_file(self):
        if self.json_file is not None:
            os.unlink(self.json_file)



class ReducedCluster(object):
    """docstring for ReducedCluster"""
    def __init__(self, seq_id):
        super(ReducedCluster, self).__init__()
        self.name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        self.sub_clusters = ['_'.join(seq_id.split('_')[:-1]), ]
        self._seq_ids = None
        self._size = None


    @property
    def size(self):
        self._size = len(self.seq_ids)
        return self._size


    @property
    def seq_ids(self):
        if self._seq_ids is None:
            return []
        return self._seq_ids

    @seq_ids.setter
    def seq_ids(self, seq_ids):
        self._seq_ids = seq_ids


    def add(self, seq_id):
        self.sub_clusters.append('_'.join(seq_id.split('_')[:-1]))


    def get_sequences(self, cluster_db, mr_db):
        seq_ids = []
        for sc in self.sub_clusters:
            seq_ids.extend(list(cluster_db.find(sc)))
        return list(mr_db.find(seq_ids))





def get_clusters_from_file(cluster_file, mr_db=None):
    '''
    From clonify_map() output (a cluster file), return a list of
    Cluster objects.
    '''
    clusters = {}
    with open(cluster_file) as f:
        for line in f:
            s, c = line.strip().split()
            clusters[c] = clusters[c] + [s] if c in clusters else [s]
    return [Cluster(seq_ids=v, mrdb=mr_db) for v in clusters.values()]



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
