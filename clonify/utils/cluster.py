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


import os
import random
import string
import tempfile


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

    @property
    def name(self):
        if self._name is None:
            self._name = ''.join(random.choice(string.ascii_uppercase + string.digits) for _ in range(8))
        return self._name


    def json(self, other=None, as_file=False, temp_dir=None):
        other = other if other is not None else []
        sequences = self.sequences + other.sequences
        if as_file:
            self.json_file = self.pretty_json(sequences, as_file, temp_dir)
            return self.json_file
        else:
            return self.pretty_json(sequences, as_file, temp_dir)


    @staticmethod
    def pretty_json(sequences, as_file=False, temp_dir=None):
        jsons = []
        for s in sequences:
            mut_list = [BASE_MUT.format(loc=m['loc'], mut=m['mut']) for m in s['var_muts_nt']['muts']]
            mut_string = ', \n'.join(mut_list)
            jsons.append(BASE_JSON.format(v_all=s['v_gene']['full'].split('*')[-1],
                                          v_gene='-'.join(s['v_gene']['full'].split('*')[0].split('-')[1:]),
                                          v_full=s['v_gene']['full'],
                                          v_fam=s['v_gene']['full'].split('-')[0].replace('IGHV', ''),
                                          seq_id=s['seq_id'],
                                          j_all=s['j_gene']['full'].split('*')[-1],
                                          j_full=s['j_gene']['full'],
                                          j_gene=s['j_gene']['full'].split('*')[0].replace('IGHJ', ''),
                                          junc_aa=s['junc_aa'],
                                          mut_string=mut_string,
                                          mut_num=s['var_muts_nt']['num']))
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
