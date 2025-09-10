from .annotations_func import isADARFixable, isApoBecFixable, vep_annotations
from .pre_process import pre_process

instructions = {
    "name": "default",
    "key_cols": ["chr", "pos", "ref", "alt"],
    "description": "Default database for variant annotations",
    "annotations": {
        # "isADARFixable": {
        #     "data_type": "bool",
        #     "description": "Indicates if the variant is ADAR fixable",
        #     "compute_function": isADARFixable,
        # },
        # "isApoBecFixable": {
        #     "data_type": "bool",
        #     "description": "Indicates if the variant is ApoBec fixable",
        #     "compute_function": isApoBecFixable,
        # },
        "vep_annotations": {
            "data_type": "list",
            "description": "VEP annotations for the variant",
            "compute_function": vep_annotations,
        },

        # Add more annotations as needed
    },
    "pre_processor": pre_process,
}