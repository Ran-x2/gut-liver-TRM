import multiprocess as mp
mp.set_start_method('fork', force=True)
from multiprocessing import freeze_support

import sys
import os
import pickle
from datetime import datetime
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
from geneformer import TranscriptomeTokenizer
from geneformer import Classifier

# f = open("250418_test_output.txt", "a")
storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/gdT_geneformer_TRMTeff'
output_prefix="gdT_activation"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"
six_digit_date = '250411'
model = f"{storage_dir}/{six_digit_date}_geneformer_cellClassifier_{output_prefix}/ksplit1/run-ca434776/checkpoint-392"


cell_states_to_model = {
    "state_key": "general_type", 
    "start_state": "TRM", 
    "goal_state": "Teff"
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

isp = InSilicoPerturber(perturb_type="overexpress",
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
                        forward_batch_size=16,
                        nproc=60)

isp.perturb_data(
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/",
    "TRM2Teff"
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
    "TRM2Teff"
)