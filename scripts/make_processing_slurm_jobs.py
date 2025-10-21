"""
Make task list and slurm scripts for processing.
"""
__author__ = "Zhe.Feng@pnnl.gov"

import yaml
import textwrap
import subprocess

def load_config(config_file, catalog_source):
    """
    Load configuration from YAML file for a specific catalog source.
    
    Args:
        config_file: str
            Path to the YAML configuration file
        catalog_source: str
            The catalog source key to load configuration for
            
    Returns:
        dict: Configuration dictionary for the specified source
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)
    
    if catalog_source not in config:
        raise ValueError(f"Catalog source '{catalog_source}' not found in config file. "
                        f"Available sources: {list(config.keys())}")
    
    return config[catalog_source]

if __name__ == "__main__":

    # Submit slurm job
    submit_job = True

    root_dir = "/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/"
    code_dir = f"{root_dir}scripts/"
    config_dir = f"{root_dir}config/"
    slurm_dir = f"{root_dir}slurm/"
    config_sources_file = f"{config_dir}config_sources.yaml"

    slurm_basename = f'slurm_'

    # Specify which code function to run
    # code_func = "make_mcs_swath_masks"
    # code_func = "combine_tracking_masks"
    # code_func = "make_cooccurrence_masks"
    code_func = "calc_monthly_rainmap_by_cof"

    # Build full code name
    code_name = f"{code_dir}{code_func}.py"

    # Set wallclock_time based on which code is used
    if "make_mcs_swath_masks" in code_name:
        wallclock_time = "00:15:00"
    elif "combine_tracking_masks" in code_name:
        wallclock_time = "00:05:00"
    elif "make_cooccurrence_masks" in code_name:
        wallclock_time = "00:15:00"
    elif "calc_monthly_rainmap_by_cof" in code_name:
        wallclock_time = "00:05:00"
    else:
        wallclock_time = "00:30:00"

    # Define a list of run names to process
    runnames = [
        "casesm2_10km_nocumulus",
        "icon_d3hp003",
        "ifs_tco3999_rcbmf",
        "IR_IMERG",
        "nicam_gl11",
        "scream_ne120",
        "um_glm_n2560_RAL3p3",
    ]

    # Create the list of job tasks needed by SLURM...
    task_filename = f"{slurm_dir}tasks_{code_func}.txt"
    task_file = open(task_filename, "w")
    ntasks = 0

    # Loop over sources
    for run in runnames:

        # Load configuration for the specified source
        config = load_config(config_sources_file, run)
        source = config.get("source_name")
        # import pdb; pdb.set_trace()

        if "make_mcs_swath_masks" in code_name:
            config_file = f"{config_dir}config_mcs_tbpf_{source}.yml"
            cmd = f"python {code_name} -c {config_file}"

        elif "combine_tracking_masks" in code_name:
            config_file = f"{config_sources_file}"
            # Special case for IMERGv7
            if source == "IMERGv7":                
                cmd = f"python {code_dir}combine_era5_imerg_tracking_masks.py"
            else:
                cmd = f"python {code_name} -c {config_file} --source {run}"

        elif "make_cooccurrence_masks" in code_name:
            config_file = f"{config_sources_file}"
            cmd = f"python {code_name} -c {config_file} --source {run}"

        elif "calc_monthly_rainmap_by_cof" in code_name:
            config_file = f"{config_sources_file}"
            cmd = f"python {code_name} -c {config_file} --source {run}"

        task_file.write(cmd + "\n")
        ntasks += 1

    task_file.close()
    print(task_filename)


    # Create a SLURM submission script for the above task list...
    slurm_filename = f'{slurm_dir}{slurm_basename}{code_func}.sh'
    slurm_file = open(slurm_filename, "w")
    text = f"""\
        #!/bin/bash
        #SBATCH -A m1867
        #SBATCH -J {code_func}
        #SBATCH -t {wallclock_time}
        #SBATCH -q regular
        #SBATCH -C cpu
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=128
        #SBATCH --exclusive
        #SBATCH --output=log_{code_func}_%A_%a.log
        #SBATCH --mail-type=END
        #SBATCH --mail-user=zhe.feng@pnnl.gov
        #SBATCH --array=1-{ntasks}

        date
        # Activate Python environment
        source activate /global/common/software/m1867/python/pyflex

        # cd {code_dir}

        # Takes a specified line ($SLURM_ARRAY_TASK_ID) from the task file
        LINE=$(sed -n "$SLURM_ARRAY_TASK_ID"p {task_filename})
        echo $LINE
        # Run the line as a command
        $LINE

        date
    """
    slurm_file.writelines(textwrap.dedent(text))
    slurm_file.close()
    print(slurm_filename)

    # Run command
    if submit_job == True:
        cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
        print(cmd)
        subprocess.run(f'{cmd}', shell=True)