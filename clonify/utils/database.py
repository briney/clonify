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


import cPickle as pickle
import os
import sqlite3


class Database(object):
    """docstring for Database"""
    def __init__(self, name, direc):
        super(Database, self).__init__()
        self.name = name
        self.path = os.path.join(direc, name)
        self.table_name = 'seqs'
        self.structure = [('key', 'text'), ('value', 'text')]
        self._connection = None
        self._cursor = None
        self._create_table_cmd = None
        self._insert_cmd = None
        self.create_table()


    @property
    def connection(self):
        if self._connection is None:
            self._connection = sqlite3.connect(self.path)
        return self._connection

    @property
    def cursor(self):
        if self._cursor is None:
            self._cursor = self.conn.cursor()
        return self._cursor

    @property
    def create_table_cmd(self):
        if self._create_cmd is None:
            field_string = ', '.join([' '.join(s) for s in structure])
            self._create_cmd = 'CREATE TABLE {} {}'.format(self.table_name,
                self.structure)
        return self._create_cmd

    @property
    def insert_cmd(self):
        if self._insert_cmd is None:
            self._insert_cmd = 'INSERT INTO {} VALUES ({})'.format(self.table_name,
                ','.join(['?'] * len(structure)))
        return self._insert_cmd


    def create_table(self):
        self.cursor.execute(self.create_table_cmd)


    def insert_one(self, value):
        '''
        Inserts a single key/value pair

        Inputs:
          value - an iterable (list or tuple) containing a single key/value pair
        '''
        if len(value) != len(self.structure):
            err = 'Mismatch between the supplied values:\n{}\n'.format(values)
            err += 'and the structure of the table:\n{}'.format(self.structure)
            raise RuntimeError(err)
        self.cursor.execute(self.insert_cmd, value)


    def insert_many(self, values):
        '''
        Inserts multiple key/value pairs.

        Inputs:
          values - a list of iterables (lists or tuples) with each iterable containing
                a single key/value pair.
        '''
        for val in values:
            if len(val) != len(self.structure):
                err = 'Mismatch between one of the supplied values:\n{}\n'.format(val)
                err += 'and the structure of the table:\n{}'.format(self.structure)
                raise RuntimeError(err)
        self.cursor.executemany(self.insert_cmd, values)


    def find(self, keys, unpickle=False):
        '''
        Searches a SQLite database.

        Inputs:
          keys - a single key (as a string) or one or more keys (as a list or tuple),
                containing the keys or which values will be returned
          unpickle - boolean, whether to unpickle (via pickle.loads()) the results
                before returning.

        Returns: a list of values
        '''
        if type(keys) in [str, unicode]:
            keys = [keys, ]
        results = []
        for chunk in self.chunker(keys):
            result_chunk = self.cursor.execute(
                '''SELECT seqs.key, seqs.value
                FROM seqs
                WHERE seqs.key IN ({})'''.format(','.join('?' * len(chunk))), chunk)
            results.extend(result_chunk)
        if unpickle:
            return [pickle.loads(str(r[1])) for r in results]
        return [r[1] for r in results]


        # # convert the projection to a SQL SELECT string
        # if project is not None:
        #     select_list = self._parse_projection(project)
        # else:
        #     select_list = ['{}.{}'.format(self.table_name, s[0]) for s in self.structure]
        # select_string = ', '.join(select_list)
        # # convert the match into a SQL WHERE string and run query
        # query = 'SELECT {} FROM {}'.format(select_string, self.table_name)
        # if match is None:
        #     results = self.cursor.execute(query)
        # else:
        #     where_list, where_vals = self._parse_match(match)
        #     where_string = ' AND '.join(where_list)
        #     query += ' WHERE {}'.format(where_string)
        #     results = self.cursor.execute(query, where_vals)
        # # fetch the results
        # return results.fetchone() if find_one else results.fetchall()


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
        Yield successive n-sized chunks from l.
        '''
        for i in xrange(0, len(l), n):
            yield l[i:i + n]


    # def _parse_projection(self, project):
    #     if type(project) in [str, unicode]:
    #         select_fields = ['{}.{}'.format(self.table_name, project)]
    #     elif type(project) in [list, tuple]:
    #         select_fields = ['{}.{}'.format(self.table_name, p) for p in project]
    #     elif type(project) == dict:
    #         select_fields = ['{}.{}'.format(self.table_name, k) for k, v in project.items() if v]
    #     else:
    #         select_fields = []
    #     return select_fields


    # def _parse_match(self, match):
    #     if type(match) != dict:
    #         err = '<match> must be a dictonary. You provided:\n{}\n'.format(match)
    #         err += 'which is {}'.format(type(match))
    #         raise TypeError(err)
    #     where_list = []
    #     where_vals = []
    #     for k, v in match.items():
    #         if type(v) in [str, unicode]:
    #             v = [v]
    #         where_vals += v
    #         where = '{} IN ({})'.format(k, ['?'] * len(v))
    #         where_list.append(where)
    #     return where_list, where_vals
