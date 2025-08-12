#srun -p gpu --gres=gpu:1 --cpus-per-task=24 --mem=128G  --time=4200 --pty /bin/bash
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

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/gdT_geneformer_TRMTeff'
output_prefix="gdT_activation"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"

# tk = TranscriptomeTokenizer({"general type": "general_type"}, nproc=39)
# tk.tokenize_data(f"{storage_dir}", 
#                  f"{storage_dir}",
#                  "tokenized", 
#                  file_format="h5ad")


# cc = Classifier(classifier="cell",
#                 cell_state_dict = {"state_key": "general_type", "states": "all"},
#                 max_ncells=None,
#                 freeze_layers = 6,
#                 num_crossval_splits = 1,
#                 split_sizes = {"train": 0.6, "valid": 0.2, "test": 0.2},
#                 forward_batch_size=16,
#                 nproc=39)

# cc.prepare_data(input_data_file=f"{storage_dir}/tokenized.dataset",
#                 output_directory=f"{storage_dir}/",
#                 output_prefix=output_prefix)

# all_metrics = cc.validate(model_directory=vanilla_model,
#                           prepared_input_data_file=f"{storage_dir}/{output_prefix}_labeled_train.dataset",
#                           id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
#                           output_directory=f"{storage_dir}/",
#                           output_prefix=output_prefix,
#                           #n_hyperopt_trials=1,
#                           predict_eval=True)

training_args = {
    "num_train_epochs": 1,
    "learning_rate": 0.0005677283501766975,
    "lr_scheduler_type": "polynomial",
    "warmup_steps": 118.9506,
    "weight_decay":0.0023215567076607796,
    "per_device_train_batch_size": 12,
    "seed": 69,
}

cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "general_type", "states": "all"},
                max_ncells=None,
                training_args=training_args,
                freeze_layers = 6,
                num_crossval_splits = 1,
                split_sizes = {"train": 0.6, "valid": 0.2, "test": 0.2},
                forward_batch_size=16,
                nproc=30)


# now = datetime.now()
# six_digit_date = now.strftime("%y%m%d")

model = f"{storage_dir}/250411_geneformer_cellClassifier_{output_prefix}/ksplit1/run-ca434776/checkpoint-392"

all_metrics = cc.evaluate_saved_model(
        model_directory=model,
        id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
        test_data_file=f"{storage_dir}/{output_prefix}_labeled_test.dataset",
        output_directory=f"{storage_dir}",
        output_prefix=output_prefix,
    )

# with open(f"{storage_dir}/{six_digit_date}_geneformer_cellClassifier_{output_prefix}/{output_prefix}_eval_metrics_dict.pkl", 'rb') as file:
#     all_metrics = pickle.load(file)

embex = EmbExtractor(model_type="CellClassifier",
                     num_classes=2, 
                     max_ncells=10000,
                     emb_layer=-1, 
                     emb_label=["general_type"],
                     labels_to_plot=["general_type"],
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

cc.plot_conf_mat(
        conf_mat_dict={"Geneformer": all_metrics["conf_matrix"]},
        output_directory=f"{storage_dir}/",
        output_prefix=output_prefix
)