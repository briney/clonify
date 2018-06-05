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
import sqlite3
import sys
import time

if sys.version_info[0] > 2:
    STR_TYPES = [str, ]
    import pickle
else:
    STR_TYPES = [str, unicode]
    import cPickle as pickle


class Database(object):
    """docstring for Database"""
    def __init__(self, name=None, direc=None, in_memory=False):
        super(Database, self).__init__()
        self.name = name
        if all([name is not None, direc is not None]):
            self.dir = os.path.abspath(direc)
            self.path = os.path.join(self.dir, self.name)
        elif in_memory:
            self.dir = None
            self.path = ':memory:'
        self.table_name = 'seqs'
        self.structure = [('key', 'text'), ('value', 'text')]
        self.initialized = False
        self._connection = None
        self._cursor = None
        self._create_table_cmd = None
        self._insert_cmd = None
        if not self.initialized:
            self.create_table()

    def __getitem__(self, key):
        return self.find_one(key)

    def __setitem__(self, key, value):
        return self.insert_one(key, value)


    @property
    def connection(self):
        if self._connection is None:
            self._connection = sqlite3.connect(self.path)
        return self._connection

    @connection.setter
    def connection(self, connection):
        self._connection = connection


    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.connection.cursor()
        return self._cursor

    @cursor.setter
    def cursor(self, cursor):
        self._cursor = cursor


    @property
    def create_table_cmd(self):
        if self._create_table_cmd is None:
            field_string = ', '.join([' '.join(s) for s in self.structure])
            self._create_table_cmd = 'CREATE TABLE {} ({})'.format(self.table_name,
                field_string)
        return self._create_table_cmd


    @property
    def insert_cmd(self):
        if self._insert_cmd is None:
            self._insert_cmd = 'INSERT INTO {} VALUES ({})'.format(self.table_name,
                ','.join(['?'] * len(self.structure)))
        return self._insert_cmd


    def commit(self):
        self.connection.commit()


    def close(self):
        self.connection.close()
        del self._cursor
        del self._connection
        self._cursor = None
        self._connection = None


    def create_table(self):
        self.cursor.execute('DROP TABLE IF EXISTS {}'.format(self.table_name))
        self.cursor.execute(self.create_table_cmd)
        self.initialized = True


    def insert_one(self, key, value):
        '''
        Inserts a single key/value pair

        Inputs:
          key - the key
          value - the value
        '''
        with self.connection as conn:
            conn.execute(self.insert_cmd, (key, pickle.dumps(value, protocol=0)))


    def insert_many(self, data):
        '''
        Inserts multiple key/value pairs.

        Inputs:
          data - a list/generator of iterables (lists or tuples) with each iterable containing
                a single key/value pair.
        '''
        # for kv in data:
        #     if len(kv) != len(self.structure):
        #         err = 'Mismatch between one of the supplied values:\n{}\n'.format(kv)
        #         err += 'and the structure of the table:\n{}'.format(self.structure)
        #         raise RuntimeError(err)
        data = ((d[0], pickle.dumps(d[1], protocol=0)) for d in data)
        with self.connection as conn:
            conn.executemany(self.insert_cmd, data)


    def find_one(self, key):
        '''
        Searches a SQLite database.

        Inputs:
          key - a single key (as a string)

        Returns: a single unpickled value
        '''
        self.cursor.execute(
            '''SELECT seqs.key, seqs.value
               FROM seqs
               WHERE seqs.key LIKE ?''', (key, ))
        return pickle.loads(self.cursor.fetchone()[1])


    def find(self, keys):
        '''
        Searches a SQLite database.

        Inputs:
          keys - a single key (string) or iterable (list/tuple) containing one or more keys

        Returns: a list of unpickled values
        '''
        if type(keys) in STR_TYPES:
            keys = [keys, ]
        results = []
        for chunk in self.chunker(keys):
            result_chunk = self.cursor.execute(
                '''SELECT seqs.key, seqs.value
                   FROM seqs
                   WHERE seqs.key IN ({})'''.format(','.join('?' * len(chunk))), chunk)
            results.extend(result_chunk)
        return [pickle.loads(r[1]) for r in results]


    def find_all(self):
        '''
        Returns all values in a SQLite database.

        Inputs:
          keys - a single key (string) or iterable (list/tuple) containing one or more keys

        Returns: a list of unpickled values
        '''
        if type(keys) in STR_TYPES:
            keys = [keys, ]
        results = self.cursor.execute(
            '''SELECT seqs.value
               FROM seqs''')
        return [pickle.loads(r[0]) for r in results]


    def delete(self, keys):
        if type(keys) in STR_TYPES:
            keys = [keys, ]
        with self.connection as conn:
            conn.executemany(
                '''DELETE
                   FROM seqs
                   WHERE seqs.key == ?''', keys)


    def index(self, field='key'):
        '''
        Indexes the database

        Inputs:
          field - the field on which to create the index. Default is 'key'.
        '''
        index_name = field + '_index'
        self.cursor.execute('CREATE INDEX {} ON {} ({})'.format(index_name,
            self.table_name,
            field))


    @staticmethod
    def chunker(l, n=900):
        '''
        Yields successive n-sized chunks from l.
        '''
        for i in range(0, len(l), n):
            yield l[i:i + n]
