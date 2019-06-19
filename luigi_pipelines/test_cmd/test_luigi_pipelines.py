import sys
from os.path import dirname

sys.path.insert(0, dirname(dirname(dirname(__file__))))
from luigi_pipelines import run_cmd
import click

project_root_path = dirname(dirname(dirname(__file__)))


@click.group()
def cli():
    " For group below function together"
    pass


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_germline(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py main_entry --tab {project_root_path}/test_set/germline/data_input.tsv --odir {odir} --analysis-type germline --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_germline_gatk4(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py main_entry --tab {project_root_path}/test_set/germline/data_input.tsv --odir {odir} --analysis-type germline_gatk4 --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_somatic(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py main_entry --tab {project_root_path}/test_set/somatic/data_input.tsv --odir {odir} --analysis-type somatic --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


@cli.command()
@click.option("-o", "--odir", help="output directory for testing ...")
def test_somatic_gatk4(odir):
    run_cmd(
        f"python3 {project_root_path}/luigi_pipelines/main.py main_entry --tab {project_root_path}/test_set/somatic/data_input.tsv --odir {odir} --analysis-type somatic_gatk4 --workers 5 --log-path {odir}/cmd_log.txt",
        dry_run=False)


if __name__ == '__main__':
    cli()
