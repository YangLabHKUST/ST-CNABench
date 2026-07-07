from .gt import GTLoader
from .infercnv_expr import InferCNVExprLoader
from .infercnv_cna import InferCNVCNALoader
from .copykat import CopyKATLoader
#from .copykat_gene import CopyKATGeneLoader
from .scevan_expr import SCEVANExprLoader
from .scevan_cna import SCEVANCNALoader
from .clonalscope_wgs import ClonalscopeWGSLoader
from .clonalscope_nowgs import ClonalscopeNoWGSLoader
from .numbat_expr import NumbatExprLoader
from .numbat_cna import NumbatCNALoader
from .calicost import CalicoSTCNALoader
from .starch import STARCHLoader
from .xclone_cna import XcloneCNALoader
from .xclone_expr import XcloneExprLoader
from .slidecna import SlideCNALoader

LOADER_REGISTRY = {
    'InferCNV_expr': InferCNVExprLoader,
    'InferCNV_cna': InferCNVCNALoader,
    'CopyKAT': CopyKATLoader,
    'SCEVAN_expr': SCEVANExprLoader,
    'SCEVAN_cna': SCEVANCNALoader,
    'Clonalscope_WGS': ClonalscopeWGSLoader,
    'Clonalscope_NoWGS': ClonalscopeNoWGSLoader,
    'Numbat_cna': NumbatCNALoader,
    'Numbat_expr': NumbatExprLoader,
    'CalicoST': CalicoSTCNALoader,
    'STARCH': STARCHLoader,
    'Xclone_cna': XcloneCNALoader,
    'Xclone_expr': XcloneExprLoader,
    'SlideCNA': SlideCNALoader,
}

EFFICIENCY_LOADER_LIST = [
    'InferCNV_expr',
    'CopyKAT',
    'SCEVAN_expr',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
    'Numbat_cna',
    'CalicoST',
    'Xclone_cna',
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
    'Numbat_cna',
    'SCEVAN_cna',
    'CopyKAT',
    'InferCNV_cna',
    'Xclone_cna',
    'STARCH',
]

# Should run in different way
TN_PRED_LOADER_LIST = [
    'SCEVAN_cna',
    'CopyKAT',
    'Numbat_cna',
    'STARCH',
    'CalicoST',
    'InferCNV_expr',
    'Xclone_cna',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
]

SUBCLONE_PRED_LOADER_LIST = [
    'InferCNV_cna',
    'CopyKAT',
    'SCEVAN_cna',
    'Clonalscope_WGS',
    'Clonalscope_NoWGS',
    'Numbat_cna',
    'CalicoST',
    'Xclone_cna',
    'STARCH',
]

SUBCLONE_IN_SLICE_LOADER_LIST = SUBCLONE_PRED_LOADER_LIST + [
    'SlideCNA',
]
