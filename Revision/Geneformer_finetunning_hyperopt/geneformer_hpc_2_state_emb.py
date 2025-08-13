import multiprocess as mp
mp.set_start_method('fork', force=True)
from multiprocessing import freeze_support

import sys
import os
import pickle
from datetime import datetime
sys.path.append(os.getcwd())
import glob
import json
from pathlib import Path
import re
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
from geneformer import TranscriptomeTokenizer
from geneformer import Classifier

# f = open("250418_test_output.txt", "a")
storage_dir = os.getcwd()
output_prefix = storage_dir.split('/mnt/vstor/SOM_PATH_DKB50/members/rxr456/Trapecar_geneformer_finetune/')[1] + '_expansion'
print(output_prefix)
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"
json_filename = glob.glob(f"{storage_dir}/*geneformer_cellClassifier*/ksplit1/_objective*/*experiment_state*.json")
print(json_filename)
text = Path(json_filename[0]).read_text(encoding="utf-8", errors="ignore")

NUM = r'[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?'
pat_acc_escaped = re.compile(r'\\"eval_macro_f1\\":\s*(' + NUM + r'),\\n')
pat_res_escaped = re.compile(
    r'\\"trial_id\\":\s*\\"(?P<val>(?:[^\\]|\\(?!"))*)\\",\\n'
)
acc_hits = [float(m) for m in pat_acc_escaped.findall(text)]
res_hits = [m for m in pat_res_escaped.findall(text)]
unique_res_hits = list(dict.fromkeys(res_hits))

max_acc_index = max(range(len(acc_hits)), key=lambda i: acc_hits[i])
trial_name = unique_res_hits[max_acc_index]

checkpoint_parent_dir = json_filename[0].split('_objective')[0]+'run-'+trial_name
subfolders = [f for f in os.listdir(checkpoint_parent_dir) if os.path.isdir(os.path.join(checkpoint_parent_dir, f))]

if subfolders:
    only_subfolder_path = os.path.join(checkpoint_parent_dir, subfolders[0])
    print(only_subfolder_path)
else:
    print("No subfolder found.")
model = only_subfolder_path

cell_states_to_model = {
    "state_key": "top10_or_not", 
    "start_state": "False", 
    "goal_state": "True"
}

# embex = EmbExtractor(model_type="CellClassifier",
#                      num_classes=2, 
#                      max_ncells=50000,
#                      emb_layer=-1, 
#                      summary_stat="exact_mean",  # I don't want this stat
#                      forward_batch_size=16,
#                      nproc=60)

# state_embs_dict = embex.get_state_embs(
#     cell_states_to_model,
#     model,
#     f"{storage_dir}/tokenized.dataset",
#     f"{storage_dir}",
#     "state_emb"
# )

with open(f"{storage_dir}/state_emb.pkl", 'rb') as file:
    state_embs_dict = pickle.load(file)

isp = InSilicoPerturber(perturb_type="overexpression",
                        genes_to_perturb="all",
                        combos=0,
                        anchor_gene=None,
                        model_type="CellClassifier",
                        num_classes=2,
                        cell_states_to_model=cell_states_to_model,
                        state_embs_dict=state_embs_dict,
                        emb_mode="cls",
                        max_ncells=5000,
                        emb_layer=-1,
                        forward_batch_size=1,
                        nproc=60)

isp.perturb_data(
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/",
    output_prefix
)

ispstats = InSilicoPerturberStats(mode="goal_state_shift",
                                  genes_perturbed="all",
                                  combos=0,
                                  anchor_gene=None,
                                  cell_states_to_model=cell_states_to_model)

ispstats.get_stats(
    f"{storage_dir}",
    None,
    f"{storage_dir}",
    output_prefix
)