"""
Make task list and slurm scripts for processing a single source with multiple code functions.
"""
__author__ = "Zhe.Feng@pnnl.gov"

import yaml
import textwrap
import subprocess
import os

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
    submit_job = False
    
    # Create bash script for interactive node execution
    make_bash_script = True

    root_dir = "/global/homes/f/feng045/program/waccem/gsrm_cooccurrencefeature/"
    code_dir = f"{root_dir}scripts/"
    config_dir = f"{root_dir}config/"
    slurm_dir = f"{root_dir}slurm/"
    log_dir = f"{slurm_dir}log/"
    config_sources_file = f"{config_dir}config_sources.yaml"

    # Create log directory if it doesn't exist
    os.makedirs(log_dir, exist_ok=True)

    slurm_basename = f'slurm_alltasks_'

    # Specify which run to process
    # runname = "scream_ne120"
    # runname = "casesm2_10km_nocumulus"
    # runname = "icon_d3hp003"
    # runname = "ifs_tco3999_rcbmf"
    runname = "IR_IMERG"
    # runname = "nicam_gl11"
    # runname = "um_glm_n2560_RAL3p3"

    # Define a list of code functions to process
    code_funcs = [
        "make_mcs_swath_masks",
        "combine_tracking_masks",
        "make_cooccurrence_masks",
        "calc_monthly_rainmap_by_cof",
    ]

    # Dictionary of wallclock times for each code function
    wallclock_times = {
        "make_mcs_swath_masks": "00:30:00",
        "combine_tracking_masks": "00:05:00",
        "make_cooccurrence_masks": "00:15:00",
        "calc_monthly_rainmap_by_cof": "00:10:00",
    }

    # Load configuration for the specified source
    config = load_config(config_sources_file, runname)
    source = config.get("source_name")

    # Create the list of job tasks needed by SLURM...
    task_filename = f"{slurm_dir}tasks_all_{runname}.txt"
    task_file = open(task_filename, "w")
    ntasks = 0

    # Loop over code functions
    for code_func in code_funcs:

        # Build full code name
        code_name = f"{code_dir}{code_func}.py"

        if code_func == "make_mcs_swath_masks":
            config_file = f"{config_dir}config_mcs_tbpf_{source}.yml"
            cmd = f"python {code_name} -c {config_file}"

        elif code_func == "combine_tracking_masks":
            config_file = f"{config_sources_file}"
            # Special case for IMERGv7
            if source == "IMERGv7":                
                cmd = f"python {code_dir}combine_era5_imerg_tracking_masks.py"
            else:
                cmd = f"python {code_name} -c {config_file} --source {runname}"

        elif code_func == "make_cooccurrence_masks":
            config_file = f"{config_sources_file}"
            cmd = f"python {code_name} -c {config_file} --source {runname}"

        elif code_func == "calc_monthly_rainmap_by_cof":
            config_file = f"{config_sources_file}"
            cmd = f"python {code_name} -c {config_file} --source {runname}"

        task_file.write(cmd + "\n")
        ntasks += 1

    task_file.close()
    print(task_filename)

    # Use the longest wallclock time from all code functions
    max_wallclock_time = max([wallclock_times.get(func, "00:30:00") for func in code_funcs])

    # Create a SLURM submission script for the above task list...
    slurm_filename = f'{slurm_dir}{slurm_basename}{runname}.sh'
    slurm_file = open(slurm_filename, "w")
    text = f"""\
        #!/bin/bash
        #SBATCH -A m1867
        #SBATCH -J {runname}
        #SBATCH -t {max_wallclock_time}
        #SBATCH -q regular
        #SBATCH -C cpu
        #SBATCH --nodes=1
        #SBATCH --ntasks-per-node=128
        #SBATCH --exclusive
        #SBATCH --output={log_dir}log_{runname}_%A_%a.log
        #SBATCH --mail-type=END
        #SBATCH --mail-user=zhe.feng@pnnl.gov
        #SBATCH --array=1-{ntasks}

        date
        # Activate Python environment
        source activate /global/common/software/m1867/python/pyflex-dev

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

    # Create bash script for interactive node execution
    if make_bash_script == True:
        bash_filename = f'{slurm_dir}run_interactive_{runname}.sh'
        bash_file = open(bash_filename, "w")
        bash_text = f"""\
            #!/bin/bash
            # Bash script to run tasks sequentially on an interactive node
            # Usage: 
            #   1. Request an interactive node:
            #      salloc -N 1 -C cpu -q interactive -t 04:00:00 -A m1867
            #   2. Run this script:
            #      bash {bash_filename}
            
            echo "=========================================="
            echo "Starting sequential task execution for {runname}"
            echo "=========================================="
            echo ""
            
            # Activate Python environment
            echo "Activating Python environment..."
            source activate /global/common/software/m1867/python/pyflex-dev
            echo ""
            
            # Run each task sequentially
            task_number=0
            while IFS= read -r line; do
                task_number=$((task_number + 1))
                echo "=========================================="
                echo "Task ${{task_number}} of {ntasks}: $line"
                echo "=========================================="
                echo "Start time: $(date)"
                echo ""
                
                # Run the command
                $line
                
                # Check if command succeeded
                if [ $? -eq 0 ]; then
                    echo ""
                    echo "✓ Task ${{task_number}} completed successfully"
                else
                    echo ""
                    echo "✗ Task ${{task_number}} FAILED with exit code $?"
                    echo "Stopping execution."
                    exit 1
                fi
                
                echo "End time: $(date)"
                echo ""
            done < {task_filename}
            
            echo "=========================================="
            echo "All {ntasks} tasks completed successfully!"
            echo "=========================================="
            echo "Final end time: $(date)"
        """
        bash_file.writelines(textwrap.dedent(bash_text))
        bash_file.close()
        
        # Make the bash script executable
        os.chmod(bash_filename, 0o755)
        
        print(f"Interactive bash script created: {bash_filename}")

    # Run command
    if submit_job == True:
        cmd = f'sbatch --array=1-{ntasks} {slurm_filename}'
        print(cmd)
        subprocess.run(f'{cmd}', shell=True)
