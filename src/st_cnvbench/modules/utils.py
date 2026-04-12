import os
import yaml
import logging
import re
from typing import Any, Dict
from pathlib import Path

def ensure_init_file(directory: str) -> bool:
    try:
        init_path = os.path.join(directory, "__init__.py")
        if not os.path.exists(init_path):
            with open(init_path, 'w'): pass
        return True
    except OSError as e:
        logging.error(f"Failed to create __init__.py in {directory}: {e}")
        return False

def load_config(config_path: str) -> Dict[str, Any]:
    """
    Loads YAML and performs recursive variable interpolation using ${VAR} syntax.
    """
    if not os.path.exists(config_path):
        raise FileNotFoundError(f"Config not found: {config_path}")

    with open(config_path, 'r') as f:
        config = yaml.safe_load(f) or {}

    context = {}
    
    # top-level scalars
    for k, v in config.items():
        if isinstance(v, (str, int, float, bool)):
            context[k] = str(v)
    if 'project_settings' in config and isinstance(config['project_settings'], dict):
        context.update({k: str(v) for k, v in config['project_settings'].items()})

    # Self-References within Context 
    for _ in range(3): 
        for k, v in context.items():
            if '${' in v:
                for src_k, src_v in context.items():
                    pattern = f"${{{src_k}}}"
                    if pattern in v:
                        v = v.replace(pattern, src_v)
                context[k] = v

    # Interpolate Config Data
    def _interpolate(node: Any) -> Any:
        if isinstance(node, dict):
            return {k: _interpolate(v) for k, v in node.items()}
        elif isinstance(node, list):
            return [_interpolate(v) for v in node]
        elif isinstance(node, str):
            if '${' not in node: return node
            # Replace all known variables
            for k, v in context.items():
                pattern = f"${{{k}}}"
                if pattern in node:
                    node = node.replace(pattern, v)
            return node
        return node

    resolved_config = _interpolate(config)
    logging.info(f"Loaded config: {config_path}")
    return resolved_config

def create_hard_link(src: Path, dst: Path):
    "create a hard link from src to dst, replacing dst if it exists"
    if dst.exists():
        dst.unlink()
    dst.hardlink_to(src)