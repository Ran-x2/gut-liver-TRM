#srun -p gpu --gres=gpu:1 --cpus-per-task=24 --mem=128G  --time=4200 --pty /bin/bash
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

cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "top10_or_not", "states": "all"},
                max_ncells=None,
                freeze_layers = 6,
                num_crossval_splits = 1,
                split_sizes = {"train": 0.6, "valid": 0.2, "test": 0.2},
                forward_batch_size=16,
                nproc=39)

# now = datetime.now()
# six_digit_date = now.strftime("%y%m%d")
all_metrics = cc.evaluate_saved_model(
        model_directory=model,
        id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
        test_data_file=f"{storage_dir}/{output_prefix}_labeled_test.dataset",
        output_directory=f"{storage_dir}",
        output_prefix=output_prefix,
    )

cc.plot_conf_mat(
        conf_mat_dict={"Geneformer": all_metrics["conf_matrix"]},
        output_directory=f"{storage_dir}/",
        output_prefix=output_prefix
)
# with open(f"{storage_dir}/{six_digit_date}_geneformer_cellClassifier_{output_prefix}/{output_prefix}_eval_metrics_dict.pkl", 'rb') as file:
#     all_metrics = pickle.load(file)

embex = EmbExtractor(model_type="CellClassifier",
                     num_classes=2,
                     emb_layer=-1, 
                     emb_label=["top10_or_not"],
                     labels_to_plot=["top10_or_not"],
                     forward_batch_size=16,
                     nproc=39)


embs = embex.extract_embs(model,
                          f"{storage_dir}/tokenized.dataset",
                          f"{storage_dir}/",
                          output_prefix + "_embeddings_labeled")

embex.plot_embs(embs=embs,
                plot_style="heatmap",
                output_directory=f"{storage_dir}/",
                output_prefix="embeddings_heatmap")