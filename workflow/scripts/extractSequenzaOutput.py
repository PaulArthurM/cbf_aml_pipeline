
import pandas





def openSequenzaConfintsCP(file):
    """Return SAMPLE_confints_CP.txt in a pandas DataFrame"""
    return pandas.read_csv(file, sep='\t', header = 0)


#'results/{token}/sequenza/{sample, [A-Za-z0-9]+}_seqz/{sample}_segments.txt'

@click.command()
@click.option('-t', '--token', default=None, help="Token.", required=True)
def main(token):
    paths = 'results/{token}/sequenza/{sample}_seqz/{sample}_segments.txt'
    sequenzaOutput = openSequenzaConfintsCP(file)


if __name__== '__main__':
    main()
