#!/usr/local/bin/python
# filename: lineage.py

#
# Copyright (c) 2019 Bryan Briney
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


from uuid import uuid4




class Lineages():
    '''
    Representation of a set of lineages, as assigned by Clonify.
    '''

    def __init__(self, seq_ids, cluster_file):
        self.seq_ids = seq_ids
        self.cluster_ids = []
        self.lineage_dict = {}
        self._process_cluster_file(cluster_file)

    
    @property
    def lineages(self):
        return list(self.lineage_dict.values())


    def _process_cluster_file(self, cluster_file):
        with open(cluster_file) as f:
            for line in f:
                if line.strip():
                    self.cluster_ids.append(line.strip())
        for cluster_id, seq_id in zip(self.cluster_ids, self.seq_ids):
            if cluster_id in self.lineage_dict:
                self.lineage_dict[cluster_id].add(seq_id)
            else:
                self.lineage_dict[cluster_id] = Lineage(seq_id)

    
        



class Lineage():
    '''
    Representation of a single lineage, as assigned by Clonify.
    '''
    def __init__(self, seq_id):
        self.seq_ids = [seq_id, ]
        self.id = uuid4()

    
    @property
    def size(self):
        return len(self.seq_ids)

    
    def add(self, seq_id):
        self.seq_ids.append(seq_id)

    
    def write(self, direc):
        filename = os.path.join(direc, self.id)
        if os.path.isfile(filename):
            self.id = uuid4()
            filename = os.path.join(direc, self.id)
        with open(filename, 'w') as f:
            f.write('\n'.join(self.seq_ids))
        return filename

