#!/usr/local/bin/python
# filename: sequence.py


###########################################################################
#
# Copyright (c) 2014 Bryan Briney.  All rights reserved.
#
# @version: 1.0.0
# @author: Bryan Briney
# @license: MIT (http://opensource.org/licenses/MIT)
#
###########################################################################



import json



class Sequence(object):
	"""docstring for Sequence"""
	def __init__(self, seq):
		super(Sequence, self).__init__()
		self.seq = seq

	def as_json(self):
		seq_id = self.seq['seq_id']
		v_full = self.seq['v_gene']['full']
		v_fam = v_full.split('-')[0]
		v_gene = '-'.join(v_full.split('*')[0].split('-')[1:])
		v_all = v_full.split('*')[1]
		j_full = self.seq['j_gene']['full']
		j_fam = j_full.split('*')[0]
		j_all = j_full.split('*')[1]
		junc_aa = self.seq['junc_aa']
		muts = self.seq['var_muts_nt']['muts']
		mut_count = self.seq['var_muts_nt']['num']
		j = OrderedDict([
			('v_gene', {'all': v_all,
						'gene': v_gene,
						'full': v_full,
						'fam': v_fam}),
			('seq_id', seq_id),
			('j_gene', {'all': v_all,
						'full': v_full,
						'gene': v_fam}),
			('junc_aa', junc_aa),
			('var_muts_nt', {'muts': muts,
							 'num': mut_count})
		])
		return json.dumps(j) + '\n'
