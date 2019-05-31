Projet Analyse CBF-AML

[TOC]



## Mercredi 29 Mai 2019

### Work

Téléchargement des fichiers BAM du dataset EGAD00001001433 sur le serveur a l'aide de l'API EGA :

-  WXS-SJCBF016_D	EGAN00001287403	**SJCBF016_D-C0DG1ACXX.5.bam.cip**	EGAF00000879120
- WXS-SJCBF016_D	EGAN00001287403	**SJCBF016_D-C0DG1ACXX.6.bam.cip**	EGAF00000879121

Téléchargement dans :

> /data1/scratch/pamesl/projet_cbf/data/bam

Tableau des identifiants :

| EGA Accession ID | EGA Object description       |
| ---------------- | ---------------------------- |
| EGAS             | EGA Study Accession ID       |
| EGAC             | EGA DAC Accession ID         |
| EGAP             | EGA Policy Accession ID      |
| EGAN             | EGA Sample Accession ID      |
| EGAR             | EGA Run Accession ID         |
| EGAX             | EGA Experiment ID            |
| EGAZ             | EGA Analysis Accession ID    |
| EGAD             | EGA Dataset Accession ID     |
| EGAB             | EGA Submission ID            |
| EGAF             | EGA File Unique Accession ID |

Installation de Miniconda [64-bit (bash installer)](https://docs.conda.io/en/latest/miniconda.html) sur le serveur. 

Installation de Samtools 1.9  dans l'environnement *samtools_env* :

> ```bash
> $ conda install -c bioconda samtools
> ```

Changer le nom du BAM :  *c42c7160-0031-4145-a440-0472cb948d7b* en *EGAR00001347178_SJCBF016_D-C0DG1ACXX.6.bam*

### Troubleshooting

1) Erreur lors du téléchargement des données de EGAF00000879121

```
EGA > download request_EGAF00000879121 path
---------------------------------------^^^^
asg.cliche.TokenException: java.lang.NumberFormatException: For input string: "path"
EGA > download request_EGAF00000879121
Files to download in this request: 1
Start Download Process: 5 (max:15) parallel threads
Iteration 1: 1 files.
Starting download: EGAR00001347179/SJCBF016_D-C0DG1ACXX.6.bam.cip  (0/1)
Failed to create directory for /data1/scratch/pamesl/projet_cbf/data/bam
```

Création d'un fichier : c42c7160-0031-4145-a440-0472cb948d7b.cip

> `
> $ samtools view EGAR00001347178_SJCBF016_D-C0DG1ACXX.5.bam 
> `
> Affiche la version BAM : OK

> `
> $ samtools view c42c7160-0031-4145-a440-0472cb948d7b 
> `
> Affiche la version BAM : OK

### To do

- [ ] Questionner sur l'emploi de la commande qsub
- [ ] Comprendre la différence entre les deux fichiers BAM de EGAN00001287403.

## Vendredi 31 Mai 2019

### Création d'un repository Git

Le projet est maintenant sous contrôle de versionnage avec Git. Création de la *branch* principale dans le dossier ~/data.

```bash
$ git status # Affiche les fichiers tracked et les modifications a commit
$ git add # Ajoute un fichier a etre tracked ou une modification a commit
```


