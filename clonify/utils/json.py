















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