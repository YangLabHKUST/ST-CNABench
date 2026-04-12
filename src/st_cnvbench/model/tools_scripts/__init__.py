from .calicost import CalicoSTModel
from .copykat import CopyKATModel
from .infercnv import InferCNVModel
from .clonalscope_nowgs import ClonalscopeNoWGSModel
from .numbat import NumbatModel  
from .clonalscope_wgs import ClonalscopeWGSModel
from .xclone import XcloneModel
from .scevan import SCEVANModel
from .starch import STARCHModel 

# Model registry
MODEL_REGISTRY = {
    "CalicoST": CalicoSTModel,
    "CopyKAT": CopyKATModel,
    "InferCNV": InferCNVModel,
    "Clonalscope_NoWGS": ClonalscopeNoWGSModel,
    "Clonalscope_WGS": ClonalscopeWGSModel,
    "Numbat": NumbatModel,
    "Xclone": XcloneModel,
    "SCEVAN": SCEVANModel,
    "STARCH": STARCHModel,
}
