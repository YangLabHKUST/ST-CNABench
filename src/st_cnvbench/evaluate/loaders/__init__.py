from .gt import GTLoader
from .infercnv_expr import InferCNVExprLoader
from .infercnv_cnv import InferCNVCNVLoader
from .copykat import CopyKATLoader
#from .copykat_gene import CopyKATGeneLoader
from .scevan_expr import SCEVANExprLoader
from .scevan_cnv import SCEVANCNVLoader
from .clonalscope_wgs import ClonalscopeWGSLoader
from .clonalscope_nowgs import ClonalscopeNoWGSLoader
from .numbat_expr import NumbatExprLoader
from .numbat_cnv import NumbatCNVLoader
from .calicost import CalicoSTCNVLoader
from .starch import STARCHLoader
from .xclone_cnv import XcloneCNVLoader
from .xclone_expr import XcloneExprLoader

LOADER_REGISTRY = {
    'InferCNV_expr': InferCNVExprLoader,
    'InferCNV_cnv': InferCNVCNVLoader,
    'CopyKAT': CopyKATLoader,
    'SCEVAN_expr': SCEVANExprLoader,
    'SCEVAN_cnv': SCEVANCNVLoader,
    'Clonalscope_WGS': ClonalscopeWGSLoader,
    'Clonalscope_NoWGS': ClonalscopeNoWGSLoader,
    'Numbat_cnv': NumbatCNVLoader,
    'Numbat_expr': NumbatExprLoader,
    'CalicoST': CalicoSTCNVLoader,
    'STARCH': STARCHLoader,
    'Xclone_cnv': XcloneCNVLoader,
    'Xclone_expr': XcloneExprLoader,
}

EFFICIENCY_LOADER_LIST = [
    'InferCNV_expr',
    'CopyKAT',
    'SCEVAN_expr',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
    'Numbat_cnv',
    'CalicoST',
    'Xclone_cnv',
    'STARCH',
]

RESOLUTION_LOADER_LIST = [
    # Gene-level (expr)
    'InferCNV_expr',
    'Xclone_expr',
    'SCEVAN_expr',
    'Numbat_expr',
    # Seg-level
    'CalicoST',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
    'Numbat_cnv',
    'SCEVAN_cnv',
    'CopyKAT',
    'InferCNV_cnv',
    'Xclone_cnv',
    'STARCH',
]

# Should run in different way
TN_PRED_LOADER_LIST = [
    'SCEVAN_cnv',
    'CopyKAT',
    'Numbat_cnv',
    'STARCH',
    'CalicoST',
    'InferCNV_expr',
    'Xclone_cnv',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
]

SUBCLONE_PRED_LOADER_LIST = [
    'InferCNV_cnv',
    'CopyKAT',
    'SCEVAN_cnv',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
    'Numbat_cnv',
    'CalicoST',
    'Xclone_cnv',
    'STARCH',
]
