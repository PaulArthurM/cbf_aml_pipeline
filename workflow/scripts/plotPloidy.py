import os
import pandas as pd
import os.path
import plotly.express as px
import click



def extract_cellularity(file, sample):
    sequenza_out = pd.read_csv(file, sep="\t")
    sequenza_out['sample'] = pd.Series(sample, index=sequenza_out.index)
    return sequenza_out


def scatter_plot_sequenza(sequenza_data, output):
    fig = px.line(sequenza_data, x="ploidy.estimate", y="cellularity", color="sample")
    fig.write_image(output)
    #fig.show()


def create_empty_df():
    column_names = ["cellularity", "ploidy.estimate", "ploidy.mean.cn", "sample"]
    return pd.DataFrame(columns = column_names)


def load_sample_sheet(sample_sheet):
    # Load sample sheet
    sample_sheet = pd.read_csv(sample_sheet, sep=";")
    SAMPLES = sample_sheet['samples']
    SAMPLE = sample_sheet.set_index("samples", drop = False)
    return SAMPLE


def create_df(SAMPLE, dir_sequenza):
    sequenza_data = create_empty_df()
    for sample in SAMPLE["samples"]:
        file = dir_sequenza + "sequenza/{sample}_seqz/{sample}_confints_CP.txt".format(sample=sample)
        if os.path.isfile(file):
            sequenza_out = extract_cellularity(file, sample)
            sequenza_data = sequenza_data.append(pd.DataFrame(data = sequenza_out), ignore_index=True)
    return sequenza_data


@click.command()
@click.option('-o', '--output', default=None, help="Output.", required=True)
@click.option('-d', '--dir-out', default=None, required=True, help="Output directory.")
@click.option('-s', '--sample-sheet', default=None, required=True, help="Sample-sheet.")
@click.option('-z', '--dir-sequenza', default=None, help="Sequenza results.", required=True)
@click.option('-n', '--dry-run', default=False, help="Dry-run.", type=bool)
def main(output, dir_out, sample_sheet, dir_sequenza, dry_run):
    SAMPLE = load_sample_sheet(sample_sheet)
    if dry_run:
        print()
    else:
        if not os.path.exists(dir_out):
            os.mkdir(dir_out)

        sequenza_data = create_df(SAMPLE, dir_sequenza)
        scatter_plot_sequenza(sequenza_data, output)


if __name__=="__main__":
    main()
