#!/usr/local/bin/python
# filename: databases.py

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
import sys
import time

from abutils.core.sequence import Sequence
from abutils.utils.database import SQLiteDatabase
from abutils.utils.utilities import nested_dict_lookup

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
    import pickle
else:
    STR_TYPES = [str, unicode]
    import cPickle as pickle



class ClonifyDB(SQLiteDatabase):
    """
    Database for storing Clonify sequence information
    """
    
    def __init__(self, name=None, direc=None, in_memory=False, table_name=None, clustering_field='vdj_nt'):
        self.clustering_field = clustering_field
        super(ClonifyDB, self).__init__(name=name, direc=direc,
                                        in_memory=in_memory, table_name=table_name)


    @property
    def structure(self):
        return [('id', 'text'), ('vgene', 'text'), ('jgene', 'text'),
                (self.clustering_field, 'text'), ('sequence', 'text')]

    
    @property
    def fields(self):
        return [['seq_id'], ['v_gene', 'gene'], ['j_gene', 'gene'], [self.clustering_field]]

    
    @property
    def insert_cmd(self):
        return 'INSERT INTO {} VALUES ({})'.format(self.table_name,
                ','.join(['?'] * len(self.structure)))
    
    
    def insert_one(self, sequence):
        '''
        Inserts a single entry.

        Args:

            data (list or dict): Either a list containing data fields (in order) or a dict
                with column/value mappings.
        '''
        data = [nested_dict_lookup(sequence, f) for f in self.fields]
        data.append(pickle.dumps(sequence))
        with self.connection as conn:
            conn.execute(self.insert_cmd, data)


    def insert_many(self, sequences):
        '''
        Inserts multiple entries.

        Args:

            data (list): Either a nested list of data fields (in order) or a list of dicts
                with column/value mappings.
        '''

        ## TODO: sequences might be a generator (huge) and we're consuming it to make the data list.
        ## would be cool if making data was also a generator (feeding executemany) or did many 
        ## execute operations without closing the transaction between them.

        data = []
        for s in sequences:
            d = [nested_dict_lookup(s, f) for f in self.fields]
            d.append(pickle.dumps(s))
            data.append(d)
        with self.connection as conn:
            conn.executemany(self.insert_cmd, data)


    def insert(self, sequences):
        if type(sequences) in [Sequence, dict]:
            sequences = [sequences, ]
        self.insert_many(sequences)




    # def get_cdr3s_for_vj_group(self, v, j):
    #     '''

    #     '''
    #     query_str = '''SELECT {0}.id, {0}.cdr3_nt, {0}.vgene, {0}.jgene
    #                    FROM {0}
    #                    WHERE {0}.vgene LIKE ? and {0}.jgene LIKE ?'''.format(self.table_name)
    #     results = self.cursor.execute(query_str, (v, j))
    #     return [(r[0], r[1]) for r in results]


    def get_seqs_for_vj_group(self, v, j, clustering_field):
        '''

        '''
        query_str = '''SELECT {0}.id, {0}.{1}, {0}.vgene, {0}.jgene
                       FROM {0}
                       WHERE {0}.vgene LIKE ? and {0}.jgene LIKE ?'''.format(self.table_name,
                                                                             clustering_field)
        results = self.cursor.execute(query_str, (v, j))
        return [(r[0], r[1]) for r in results]

    
    def get_all_sequences(self):
        '''

        '''
        query = '''SELECT {0}.sequence
                   FROM {0}'''.format(self.table_name)
        results = self.cursor.execute(query)
        return [pickle.loads(r[0]) for r in results]


    def get_sequences_by_id(self, ids):
        '''

        '''
        if type(ids) in STR_TYPES:
            ids = [ids, ]
        results = []
        for chunk in self.chunker(ids):
            result_chunk = self.cursor.execute(
                '''SELECT {0}.id, {0}.sequence
                   FROM {0}
                   WHERE {0}.id IN ({1})'''.format(self.table_name, ','.join('?' * len(chunk))), chunk)
            results.extend(result_chunk)
        return [pickle.loads(r[1]) for r in results]


    def find_distinct(self, field):
        '''

        '''
        query_str = '''SELECT DISTINCT {0}.{1}
                       FROM {0}'''.format(self.table_name, field)
        results = self.cursor.execute(query_str)
        return [r[0] for r in results]



    def find_one(self, where):
        '''
        Searches a SQLite database.

        Args:

          key - a single key (as a string)

        Returns: 
            
            A single unpickled entry
        '''
        where_vals = []
        wheres = []
        for k, v in where.items():
            if type(v) in STR_TYPES:
                wheres.append('{}.{} LIKE ?'.format(self.table_name, k))
                where_vals.append(v)
            else:
                wheres.append('{}.{} IN ({})'.format(self.table_name, k, ','.join('?' * len(v))))
                where_vals += v
        where_str = 'WHERE ' + ' AND '.join(wheres)
        
        query_str = '''SELECT *
                       FROM {}
                       {}'''.format(select_str, self.table_name, where_str)
        self.cursor.execute(query_str, wheres)
        ret = self.cursor.fetchone()
        vals = [ret[0]] + [pickle.loads(r) for r in ret[1:]]
        return {k: v for k, v in zip(self.structure)}




    def index(self, fields='id'):
        '''
        Indexes the database

        Args:
          
            fields (str or iterable): the field or list of field on which to create the index.
                Default is ``'id'``.
        '''
        if type(fields) in STR_TYPES:
            fields = [fields]
        index_name = '__'.join(fields) + '__index'
        self.cursor.execute('CREATE INDEX {} ON {} ({})'.format(index_name,
            self.table_name,
            ', '.join(fields)))


