import os
import pandas as pd
import argparse


def calculate_score_checkV(row):
    try:
        if row['provirus'] == 'Yes':
            score1 = -1
        else:
            score1 = 0

        if float(row['viral_genes']) == 0 and float(row['host_genes']) >= 1:
            score2 = -1
        else:
            score2 = 0

        if float(row['host_genes']) > 50 and row['provirus'] == 'No':
            score3 = -1
        else:
            score3 = 0

        total_scores = [score1, score2, score3]

    except KeyError:
        total_scores = [0, 0, 0]

    return total_scores


def calculate_score_virsorter(row):
    score = -1 if float(row['length']) > 50000 and float(row['hallmark']) <= 1 else 0
    return score




def process_checkV_results(input_file, output_file):
    df_checkV = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'Some Other Column1', 'provirus', 'Some Other Column2', 'Some Other Column3', 'viral_genes', 'host_genes', 'Some Other Column4', 'Some Other Column5', 'Some Other Column6', 'Some Other Column7', 'Some Other Column8', 'Some Other Column9', 'Some Other Column10'], skiprows=1)
    scores_list = df_checkV.apply(calculate_score_checkV, axis=1)
    df_checkV['Score1'], df_checkV['Score2'], df_checkV['Score3'] = zip(*scores_list)
    df_checkV.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2', 'Score3'])


def process_virsorter2_results(input_file, output_file):
    df_virsorter2 = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'Some Other Column1', 'Some Other Column2', 'Some Other Column4', 'Some Other Column3', 'length', 'hallmark', 'Some Other Column5', 'Some Other Column6'], skiprows=1)
    df_virsorter2['Sequence ID'] = df_virsorter2['Sequence ID'].apply(lambda x: x.split('||')[0])
    df_virsorter2['Score'] = df_virsorter2.apply(calculate_score_virsorter, axis=1)
    df_virsorter2.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score'])





def calculate_score_genomad(row):
    try:
        if float(row['n_plasmid_hallmarks']) >= 1:
            score1 = -1
        else:
            score1 = 0

        if float(row['marker_enrichment_c']) > 0:
            score2 = -1
        else:
            score2 = 0

        if float(row['marker_enrichment_p']) > 0:
            score3 = -1
        else:
            score3 = 0

        total_scores = [score1, score2, score3]

    except KeyError:
        total_scores = [0, 0, 0]

    return total_scores


def process_genomad_results(input_file, output_file):
    df_geNomad = pd.read_csv(input_file, sep='\t', names=['Sequence ID', 'OColumn1', 'OColumn2', 'n_plasmid_hallmarks', 'n_virus_hallmarks', 'OColumn4', 'OColumn5', 'OColumn6', 'OColumn7', 'OColumn8', 'OColumn9', 'OColumn10', 'OColumn11', 'OColumn12', 'OColumn13', 'OColumn14', 'OColumn15', 'OColumn16', 'OColumn17', 'OColumn18', 'OColumn19', 'OColumn20', 'OColumn21', 'OColumn22', 'OColumn23', 'OColumn24', 'OColumn25', 'OColumn26', 'OColumn27', 'OColumn28', 'OColumn29', 'marker_enrichment_c', 'marker_enrichment_p','marker_enrichment_v'], skiprows=1)
    scores_list = df_geNomad.apply(calculate_score_genomad, axis=1)
    df_geNomad['Score1'], df_geNomad['Score2'], df_geNomad['Score3'] = zip(*scores_list)
    df_geNomad.to_csv(output_file, index=False, sep='\t', columns=['Sequence ID', 'Score1', 'Score2', 'Score3'])






if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Process CheckV and VirSorter2 results and calculate scores for multiple files.')
    parser.add_argument('-c', '--checkV_input_dir', type=str, help='Path to the directory containing CheckV results files')
    parser.add_argument('-v', '--virsorter2_input_dir', type=str, help='Path to the directory containing VirSorter2 results files')
    parser.add_argument('-g', '--genomad_input_dir', type=str, help='Path to the directory containing genomad results files')

    args = parser.parse_args()

    if args.checkV_input_dir:
        checkV_files = [f for f in os.listdir(args.checkV_input_dir) if f.endswith('quality_summary.tsv')]
        for checkV_file in checkV_files:
            input_file = os.path.join(args.checkV_input_dir, checkV_file)
            output_file = os.path.join(args.checkV_input_dir, os.path.splitext(checkV_file)[0] + '_checkV_rmv.tsv')  
            process_checkV_results(input_file, output_file)

    if args.virsorter2_input_dir:
        virsorter2_files = [f for f in os.listdir(args.virsorter2_input_dir) if f.endswith('final-viral-score.tsv')]
        for virsorter2_file in virsorter2_files:
            input_file = os.path.join(args.virsorter2_input_dir, virsorter2_file)
            output_file = os.path.join(args.virsorter2_input_dir, os.path.splitext(virsorter2_file)[0] + '_vs2_rmv.tsv')  
            process_virsorter2_results(input_file, output_file)

    if args.genomad_input_dir:
        genomad_files = [f for f in os.listdir(args.genomad_input_dir) if f.endswith('contigs_features.tsv')]
        for genomad_file in genomad_files:
            input_file = os.path.join(args.genomad_input_dir, genomad_file)
            output_file = os.path.join(args.genomad_input_dir, os.path.splitext(genomad_file)[0] + '_genomad_rmv.tsv')  
            process_genomad_results(input_file, output_file)
