from abc import ABC, abstractmethod
from pathlib import Path
from typing import Dict, List, Any, Optional
import subprocess
import os
import time
import shutil
import logging

class BaseModel(ABC):
    """
    Abstract Base Class (ABC) for all models.
    """

    def __init__(
        self, 
        project_root: str, 
        model_cfg: Dict[str, Any],
        result_dir: str,
        model_name: str,
        exec_mode: str,
    ):
        """
        Initialize the model wrapper.
        
        Args:
            project_root: Root directory of the project.
            model_cfg: The COMPLETE parsed model configuration dictionary.
            result_dir: Directory to store results.
            model_name: The key used to look up specific settings in model_cfg.
            exec_mode: Execution mode, either "docker", "conda", or "apptainer".
        """
        self.project_root = Path(project_root)
        self.result_root = Path(result_dir)
        self.name = model_name
        self.global_config = model_cfg
        if self.name in model_cfg:
            self.model_cfg = self.global_config[self.name]
        else:
            logging.warning(f"[{self.name}] Config section not found in global config.")
            self.model_cfg = {}
        if exec_mode == "conda":
            self.env_name = self.model_cfg.get("env_name") 
        elif exec_mode == "docker":
            self.docker_image = self.model_cfg.get("docker_image")
        elif exec_mode == "apptainer":
            self.apptainer_sif = Path(self.model_cfg.get("apptainer_sif", "")).resolve()
        self.exec_mode = exec_mode
        self.enabled = self.model_cfg.get("enabled", True)


    def is_enabled_for_dataset(self, dataset_cfg: Dict[str, Any]) -> bool:
        if not self.enabled:
            return False
            
        ds_id = dataset_cfg.get("dataset_id")
        per_ds = self.model_cfg.get("per_dataset", {})
        
        if ds_id in per_ds:
            return per_ds[ds_id].get("enabled", True)
        return True
    
    def run(self, dataset_cfg: Dict[str, Any], overwrite: bool = False, verbose:bool = False):
        """
        Main execution workflow for the model on a given dataset.
        args:
            dataset_cfg: Configuration dictionary for the dataset.
            overwrite: Whether to overwrite existing results. Defaults to False.
        """
        ds_id = dataset_cfg.get('dataset_id', 'unknown_dataset')

        if not self.is_enabled_for_dataset(dataset_cfg):
            logging.warning(f"[{self.name}] Skipped: Not enabled for dataset {ds_id}")
            return
        
        logging.info(f"[{self.name}] Start processing: {ds_id}")
        
        # output dictionary
        output_dir = self.result_root / ds_id / self.name

        # Overwrite
        if output_dir.exists():
            if overwrite:
                logging.warning(f"[{self.name}] [Overwrite] Removing existing directory: {output_dir}")
                shutil.rmtree(output_dir)  
                output_dir.mkdir(parents=True, exist_ok=True)
            else:
                logging.info(f"[{self.name}] Output directory exists. Checking if we can skip...")
        else:
            output_dir.mkdir(parents=True, exist_ok=True)
        logging.info(f"[{self.name}] Output Dir: {output_dir}")

        # time
        start_time = time.time()

        logging.info(f"[{self.name}] Step 1: Preparing inputs...")
        input_files = self.prepare_inputs(dataset_cfg, output_dir, overwrite=overwrite)

        # Skip execution
        if input_files.get("_skip_execution", False):
            logging.info(f"[{self.name}] Step 2 & 3: Skipped (Files exist and overwrite=False).")
        else:
            logging.info(f"[{self.name}] Step 2: Building command...")
            cmd = self.build_command(dataset_cfg, input_files, output_dir)

            logging.info(f"[{self.name}] Step 3: Executing...")

            run_log = output_dir / f"{self.name}_run.log"
            self._execute_in_env(cmd, cwd=self.project_root, log_file=run_log, verbose=verbose)

        duration = time.time() - start_time
        logging.info(f"[{self.name}] Finished successfully ({duration:.2f}s)")


    # Abstract Interfaces 
    @abstractmethod
    def prepare_inputs(
        self, 
        dataset_cfg: Dict[str, Any], 
        output_dir: Path,
        overwrite: bool = False
    ) -> Dict[str, Path]:
        """
        [Abstract] Convert the input dataset into model-specific physical files 
        """
        pass

    @abstractmethod
    def build_command(
        self, 
        dataset_cfg: Dict[str, Any], 
        input_files: Dict[str, Path], 
        output_dir: Path
    ) -> List[str]:
        """
        [Abstract] Construct the actual CLI command list for the model script.
        """
        pass

    def _execute_in_env(self, cmd: List[str], cwd, log_file: Path, verbose: bool = False):
        project_root_abs = str(self.project_root.resolve())
        
        # Ensure .snakemake/log directory exists with proper permissions for Docker/Apptainer
        snakemake_log_dir = self.project_root / ".snakemake" / "log"
        if not snakemake_log_dir.exists():
            snakemake_log_dir.mkdir(parents=True, exist_ok=True)
            # Set permissions to allow write access
            os.chmod(snakemake_log_dir, 0o777)
        
        # Also ensure output directory and subdirectories have proper permissions for Docker
        if self.exec_mode in ["docker", "apptainer"]:
            # Get output directory from log_file path
            output_parent = log_file.parent
            # Set permissions on output directory recursively
            for dir_path in [output_parent] + list(output_parent.rglob('*')):
                if dir_path.is_dir():
                    try:
                        os.chmod(dir_path, 0o777)
                    except Exception:
                        pass  # Ignore errors
        
        # Apptainer mode
        if self.exec_mode == "apptainer":
            if not self.apptainer_sif.exists():
                raise ValueError(f"[{self.name}] Apptainer SIF not found in config")
            exec_prefix = [
                "apptainer", "exec",
                "--cleanenv",
                "--env", "HOME=/tmp",        
                "-B", f"{project_root_abs}:{project_root_abs}", 
                "--pwd", project_root_abs,                     
                str(self.apptainer_sif)
            ]
            final_cmd = exec_prefix + cmd

        # Docker mode
        elif self.exec_mode == "docker":
            if not self.docker_image:
                raise ValueError(f"[{self.name}] Docker image not specified in config.")
            container_cmd = [str(arg) for arg in cmd]

            # docker cmd
            exec_prefix = [
                "docker", "run", "--rm",
                "--user", f"{os.getuid()}:{os.getgid()}",
                "-v", f"{project_root_abs}:{project_root_abs}", 
                "-w", project_root_abs,
                "-e", "HOME=/tmp",
                self.docker_image
            ]
            final_cmd = exec_prefix + container_cmd
        # Conda mode         
        elif self.exec_mode == "conda":
            perf_file = log_file.with_suffix(".perf")
            exec_prefix = []
            if self.env_name:
                exec_prefix = ["conda", "run", "-n", self.env_name, "--no-capture-output"]
            full_run_cmd = exec_prefix + cmd
            final_cmd = ["/usr/bin/time", "-v", "-o", str(perf_file)] + full_run_cmd

        logging.info(f"[{self.name}] Final Command: {' '.join(final_cmd)}")

        logging.info(f"[{self.name}] Log file: {log_file}")
        with open(log_file, "w") as f_log:
            process = subprocess.Popen(
                final_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                cwd=cwd,
            )

            for line in process.stdout:
                f_log.write(line) 
                if verbose:       
                    print(f"    | {line.strip()}")

            return_code = process.wait()

        if return_code != 0:
            error_tail = self._tail_log(log_file, lines=20)
            banner = "=" * 50
            logging.error(f"\n{banner}\n"
                          f"CRITICAL ERROR in {self.name} (Exit Code: {return_code})\n"
                          f"Log File: {log_file}\n"
                          f"Last 20 lines:\n"
                          f"{'-'*50}\n"
                          f"{error_tail}\n"
                          f"{banner}")
            raise RuntimeError(f"{self.name} execution failed.")
        
    def _tail_log(self, log_path: Path, lines: int = 20) -> str:
        with open(log_path, "r") as f:
            return "".join(f.readlines()[-lines:])
        
    def get_script_path(self, script_name: str, subfolder: str = None) -> Path:
        base_dir = Path(__file__).parent
        
        target = base_dir / "tools_scripts"
        if subfolder:
            target = target / subfolder
        target = target / script_name
        
        if not target.exists():
            raise FileNotFoundError(f"Script not found: {target}")
            
        return target
        