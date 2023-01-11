from enum import Enum

CL_ROOT_ID = 'CL:0000000'


class COLUMNS(Enum):
    CELL_ID = "cell_id"
    ANNOTATION = "annotation"
    PREDICTION_SCORE = "prediction_score"
    DATASET_ID = "dataset_id"
    GLOBUS_ID = "globus_id"
