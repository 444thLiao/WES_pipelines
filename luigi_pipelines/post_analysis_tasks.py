
from luigi_pipelines import run_cmd,config

project_root_path = config.project_root_path


def germline_filter(indir, odir, tab):
    run_cmd(f"python3 {project_root_path}/api/var_filters.py -i {indir} --tab {tab} -o {odir}",
            dry_run=False)


def somatic_filter():
    # if needed??
    # maybe not because of pcgr....
    pass
