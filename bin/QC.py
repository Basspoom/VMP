import os
import subprocess
import pandas as pd
import json
import argparse
import yaml
import shlex
import shutil
from rich.console import Console
from rich.panel import Panel
from rich.table import Table
from rich.rule import Rule
from rich.text import Text
from rich.markdown import Markdown


def Getfastq_list(raw_data):
    path = raw_data
    fastq_list = []
    for item in os.listdir(path):
        item_path = os.path.join(path, item)
        if os.path.isdir(item_path):
            fastq_list.append(item)
    print('Your sample list:', fastq_list)
    return fastq_list

def fastp(out_dir, raw_data, fastq_list, qualified_quality, unqualified_percent, length_required, length_limit, complexity_threshold, adapter, correction, trim_front1, trim_tail1, max_len1, trim_front2, trim_tail2, max_len2, thread):
    print("Starting quality control process...")
    qc_dir = os.path.join(out_dir, 'QC')
    os.makedirs(qc_dir, exist_ok=True)

    for current_fasta in fastq_list:
        sample_dir = os.path.join(qc_dir, current_fasta)
        os.makedirs(sample_dir, exist_ok=True)
        print(f"Processing sample: {current_fasta}")

        fastp_command = f"fastp -q {qualified_quality} -u {unqualified_percent} --length_required {length_required} --length_limit {length_limit} -Y {complexity_threshold} --unpaired1 {os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R1.fastq.gz')} --unpaired2 {os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R2.fastq.gz')}"
        if adapter:
            fastp_command += f" -a {adapter}"
        if correction is not None and correction == 'True':
            fastp_command += f" --correction"
        if trim_front1:
            fastp_command += f" -f {trim_front1}"
        if trim_tail1:
            fastp_command += f" -t {trim_tail1}"
        if max_len1:
            fastp_command += f" -b {max_len1}"
        if trim_front2:
            fastp_command += f" -F {trim_front2}"
        if trim_tail2:
            fastp_command += f" -T {trim_tail2}"
        if max_len2:
            fastp_command += f" -B {max_len2}"

        fastp_command += f" -i {os.path.join(raw_data, current_fasta, f'*_1.fastq')} -o {os.path.join(sample_dir, f'{current_fasta}_trim_R1.fastq.gz')}"
        fastp_command += f" -I {os.path.join(raw_data, current_fasta, f'*_2.fastq')} -O {os.path.join(sample_dir, f'{current_fasta}_trim_R2.fastq.gz')}"
        fastp_command += f" -h {os.path.join(sample_dir, f'{current_fasta}_fastp.html')} -j {os.path.join(sample_dir, f'{current_fasta}_fastp.json')}"

        os.system(fastp_command)

        stat_unpaired_reads(os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R1.fastq.gz'), current_fasta, sample_dir)
        stat_unpaired_reads(os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R2.fastq.gz'), current_fasta, sample_dir)
        generate_fastp_summary(os.path.join(sample_dir, f'{current_fasta}_fastp.json'), current_fasta, sample_dir)

        print(f"Sample {current_fasta} has completed quality control.")

def stat_unpaired_reads(unpaired_file, sample, output_dir):
    if os.path.exists(unpaired_file):
        result = subprocess.run(['seqkit', 'stats', '-T', unpaired_file], capture_output=True, text=True)
        stats_df = pd.DataFrame([line.split('\t') for line in result.stdout.strip().split('\n')])
        stats_df.columns = stats_df.iloc[0]  
        stats_df = stats_df[1:]  
        stats_df.insert(0, "SAMPLE", sample)  
        stats_df.to_csv(os.path.join(output_dir, f"{sample}_unpaired_summary.tsv"), sep='\t', index=False)
        print(f"Unpaired reads statistics saved for sample: {sample}")
    else:
        print(f"Unpaired file not found: {unpaired_file}")

def generate_fastp_summary(fastp_json, sample, output_dir):
    if os.path.exists(fastp_json):
        with open(fastp_json) as f:
            data = json.load(f)

        before_df = pd.DataFrame(data['summary']['before_filtering'], index=[0])
        before_df.insert(0, "SAMPLE", sample)
        before_df.to_csv(os.path.join(output_dir, f"{sample}_fastp_summary_before.tsv"), sep='\t', index=False)

        after_df = pd.DataFrame(data['summary']['after_filtering'], index=[0])
        after_df.insert(0, "SAMPLE", sample)
        after_df.to_csv(os.path.join(output_dir, f"{sample}_fastp_summary_after.tsv"), sep='\t', index=False)

        print(f"fastp summary statistics saved for sample: {sample}")
    else:
        print(f"fastp JSON file not found: {fastp_json}")



def rm_host(out_dir, fastq_list, genome_fasta, genome_index, thread, raw_data):
    print("Starting remove host reads process...")
    rm_host_dir = os.path.join(out_dir, 'rm_host')
    os.makedirs(rm_host_dir, exist_ok=True)
    os.system(f"mkdir {out_dir}/Clean_reads")

    qc_dir = os.path.join(out_dir, 'QC')

    origin_quality_dir = os.path.join(out_dir, 'Clean_reads', "origin_quality")
    filtered_quality_dir = os.path.join(out_dir, 'Clean_reads', "filtered_quality")
    os.makedirs(origin_quality_dir, exist_ok=True)
    os.makedirs(filtered_quality_dir, exist_ok=True)

    for current_fasta in fastq_list:
        sample_dir = os.path.join(rm_host_dir, current_fasta)
        os.makedirs(sample_dir, exist_ok=True)
        sample_dir2 = os.path.join(qc_dir, current_fasta)

        print(f"Processing sample: {current_fasta}")

        index_files = [
            f"{genome_index}.1.bt2",
            f"{genome_index}.2.bt2",
            f"{genome_index}.3.bt2",
            f"{genome_index}.4.bt2",
            f"{genome_index}.rev.1.bt2",
            f"{genome_index}.rev.2.bt2",
        ]

        if all(os.path.exists(index_file) for index_file in index_files):
            print("Bowtie2 index files found, skipping index building.")
        else:
            bowtie2_build = f"bowtie2-build {genome_fasta} {genome_index} --threads {thread} -q"
            try:
                subprocess.run(bowtie2_build, shell=True, check=True)
                print("Bowtie2 index built successfully.")
            except subprocess.CalledProcessError as e:
                print(f"An error occurred while building Bowtie2 index: {e}")
                continue  

        bowtie2_command = f"bowtie2 -x {genome_index} -p {thread} -1 {os.path.join(sample_dir2, f'{current_fasta}_trim_R1.fastq.gz')} -2 {os.path.join(sample_dir2, f'{current_fasta}_trim_R2.fastq.gz')} -S {os.path.join(sample_dir, f'{current_fasta}_host_mapped.sam')}"

        samtools_command1 = f"samtools view -@ {thread} -bSh {os.path.join(sample_dir, f'{current_fasta}_host_mapped.sam')} > {os.path.join(sample_dir, f'{current_fasta}_host_mapped.bam')}"
        samtools_command2 = f"samtools sort -@ {thread}  {os.path.join(sample_dir, f'{current_fasta}_host_mapped.bam')}  -o {os.path.join(sample_dir, f'{current_fasta}_host_mapped.sorted.bam')}"
        samtools_command3 = f"samtools view -@ {thread}  -b  -f 12  -F 256  {os.path.join(sample_dir, f'{current_fasta}_host_mapped.sorted.bam')} > {os.path.join(sample_dir, f'{current_fasta}_host_unmapped.bam')}"
        samtools_command4 = f"samtools sort  -O BAM  -@ {thread} -n {os.path.join(sample_dir, f'{current_fasta}_host_unmapped.bam')} -o {os.path.join(sample_dir, f'{current_fasta}_host_unmapped_sorted.bam')}"
        samtools_command5 = f"samtools fastq -@ {thread}  {os.path.join(sample_dir, f'{current_fasta}_host_unmapped_sorted.bam')}  -1 {os.path.join(sample_dir, f'{current_fasta}_host_removed_1.fastq')}  -2 {os.path.join(sample_dir, f'{current_fasta}_host_removed_2.fastq')}  -n"

        try:
            subprocess.run(bowtie2_command, shell=True, check=True)
            subprocess.run(samtools_command1, shell=True, check=True)
            subprocess.run(samtools_command2, shell=True, check=True)
            subprocess.run(samtools_command3, shell=True, check=True)
            subprocess.run(samtools_command4, shell=True, check=True)
            subprocess.run(samtools_command5, shell=True, check=True)
            print(f"Sample {current_fasta} has completed host removal process.")

            os.system(f"mkdir {out_dir}/Clean_reads/{current_fasta}")
            os.system(f"cp {os.path.join(sample_dir, f'{current_fasta}_host_removed_1.fastq')} {out_dir}/Clean_reads/{current_fasta}")
            os.system(f"cp {os.path.join(sample_dir, f'{current_fasta}_host_removed_2.fastq')} {out_dir}/Clean_reads/{current_fasta}")


            stat_unpaired_reads(os.path.join(sample_dir, f'{current_fasta}_host_unmapped_sorted.bam'), current_fasta, sample_dir)
        except subprocess.CalledProcessError as e:
            print(f"An error occurred while processing sample {current_fasta}: {e}")

        origin_sample_dir = os.path.join(origin_quality_dir, current_fasta)
        filtered_sample_dir = os.path.join(filtered_quality_dir, current_fasta)
        os.makedirs(origin_sample_dir, exist_ok=True)
        os.makedirs(filtered_sample_dir, exist_ok=True)

        fastqc_filtered_command = f"fastqc -t {thread} {out_dir}/Clean_reads/*/*.fastq -o {filtered_sample_dir}"
        fastqc_origin_command = f"fastqc -t {thread} {os.path.join(raw_data, current_fasta, '*.fastq')} -o {origin_sample_dir}"
        os.system(fastqc_filtered_command)
        os.system(fastqc_origin_command)

            

def stat_unpaired_reads(unpaired_file, sample, output_dir):
    if os.path.exists(unpaired_file):
        result = subprocess.run(['seqkit', 'stats', '-T', unpaired_file], capture_output=True, text=True)
        stats_df = pd.DataFrame([line.split('\t') for line in result.stdout.strip().split('\n')])
        stats_df.columns = stats_df.iloc[0]  
        stats_df = stats_df[1:]  
        stats_df.insert(0, "SAMPLE", sample) 
        stats_df.to_csv(os.path.join(output_dir, f"{sample}_unpaired_summary.tsv"), sep='\t', index=False)
        print(f"Unpaired reads statistics saved for sample: {sample}")
    else:
        print(f"Unpaired file not found: {unpaired_file}")



def Nonpareil(out_dir, fastq_list, thread, min_reads_cut, nonpareil_mth, max_reads, kmer_length, num_subsamples, min_overlap, max_ram, verbose, random_seed):
    print("Starting redundancy and distribution assessment process...")
    nonpareil_dir = os.path.join(out_dir, 'nonpareil')
    os.makedirs(nonpareil_dir, exist_ok=True)
    print(f"Creating output directory: {nonpareil_dir}")

    rm_host_dir = os.path.join(out_dir, 'rm_host')

    for current_fasta in fastq_list:
        sample_dir = os.path.join(nonpareil_dir, current_fasta)
        os.makedirs(sample_dir, exist_ok=True)
        sample_dir2 = os.path.join(rm_host_dir, current_fasta)
        
        fastq_1 = os.path.join(sample_dir2, f'{current_fasta}_host_removed_1.fastq')
        fastq_2 = os.path.join(sample_dir2, f'{current_fasta}_host_removed_2.fastq')
        output_fq = os.path.join(sample_dir, f'{current_fasta}_rmhost.fastq')

        if not os.path.exists(fastq_1) or not os.path.exists(fastq_2):
            print(f"FASTQ files do not exist for sample: {current_fasta}")
            continue

        cmd1 = f"seqkit seq {fastq_1} -j {thread} -m {min_reads_cut} >> {output_fq}"
        result1 = subprocess.run(cmd1, shell=True)
        if result1.returncode != 0:
            print(f"Error occurred while processing {fastq_1}")
            continue 

        cmd2 = f"seqkit seq {fastq_1} -j {thread} -m {min_reads_cut} >> {output_fq}"
        result2 = subprocess.run(cmd2, shell=True)
        if result2.returncode != 0:
            print(f"Error occurred while processing {fastq_2}")
            continue 

        nonpareil_cmd = f"nonpareil -T {nonpareil_mth} -f fastq -t {thread} -s {output_fq} -b {os.path.join(sample_dir, current_fasta)}"

        if max_reads:
            nonpareil_cmd += f" -X {max_reads}"
        if kmer_length:
            nonpareil_cmd += f" -k {kmer_length}"
        if num_subsamples:
            nonpareil_cmd += f" -n {num_subsamples}"
        if min_overlap:
            nonpareil_cmd += f" -L {min_overlap}"
        if max_ram:
            nonpareil_cmd += f" -R {max_ram}"
        if verbose:
            nonpareil_cmd += f" -v {verbose}"
        if random_seed:
            nonpareil_cmd += f" -r {random_seed}"

        result_nonpareil = subprocess.run(nonpareil_cmd, shell=True)
        if result_nonpareil.returncode != 0:
            print(f"Error occurred while running Nonpareil for {current_fasta}")
            continue

        npo_file = os.path.join(sample_dir, f'{current_fasta}.npo')
        diversity_index_file = os.path.join(sample_dir, f'{current_fasta}_nonpareil_index.tsv')

        npo_file_quoted = shlex.quote(npo_file)
        diversity_index_file_quoted = shlex.quote(diversity_index_file)

        r_script = f"""
        R -e "library(Nonpareil);
            nps <- Nonpareil.set('{npo_file_quoted}');
            npsSummary <- summary(nps);
            npsSummaryIndex <- cbind(SAMPLE = rownames(npsSummary), npsSummary);
            write.table(npsSummaryIndex, file='{diversity_index_file_quoted}', quote=FALSE, sep='\\t', row.names=FALSE);"
        """

        try:
            subprocess.run(r_script, shell=True, check=True)
            print(f"Completed processing for sample: {current_fasta}")


            pdf_file = 'Rplots.pdf' 
            current_directory = os.getcwd()  
            pdf_source_path = os.path.join(current_directory, pdf_file)

            new_pdf_file = f"{current_fasta}_nonpareil_curves.pdf"
            pdf_destination_path = os.path.join(sample_dir, new_pdf_file)

            shutil.move(pdf_source_path, pdf_destination_path)

        except subprocess.CalledProcessError as e:
            print(f"An error occurred: {e}")
        except FileNotFoundError:
            print(f"PDF file not found: {pdf_source_path}")
        except Exception as e:
            print(f"An error occurred while moving the PDF: {e}")


def KMC(out_dir, fastq_list, thread, kmer_l, min_kmer_coverage, max_kmer_coverage, kmer_times, memory, min_freq, max_freq):
    print("Starting the process of k-mer frequency statistics and genome feature assessment....")
    kmc_dir = os.path.join(out_dir, 'KMC')
    os.makedirs(kmc_dir, exist_ok=True)
    print(f"Creating output directory: {kmc_dir}")

    rm_host_dir = os.path.join(out_dir, 'rm_host')

    genomescope_dir = (CONFIG.get('genomescope2.0') if 'CONFIG' in globals() else None)
    if genomescope_dir:
        genomescope_path = os.path.join(str(genomescope_dir).strip().strip('\'\"'), 'genomescope.R')
    else:
        genomescope_path = subprocess.check_output("whereis genomescope.R | awk '{print $2}'", shell=True).decode().strip()

    for current_fasta in fastq_list:
        sample_dir = os.path.join(kmc_dir, current_fasta)
        os.makedirs(sample_dir, exist_ok=True)
        sample_dir2 = os.path.join(rm_host_dir, current_fasta)

        fastq_1 = os.path.join(sample_dir2, f'{current_fasta}_host_removed_1.fastq')
        fastq_2 = os.path.join(sample_dir2, f'{current_fasta}_host_removed_2.fastq')

        merged_fastq = os.path.join(sample_dir, f'{current_fasta}_merge.assembled.fastq')
        with open(merged_fastq, 'wb') as f_out:
            with open(fastq_1, 'rb') as f_in:
                f_out.write(f_in.read())
            with open(fastq_2, 'rb') as f_in:
                f_out.write(f_in.read())

        output_dir = os.path.join(sample_dir, current_fasta)
        os.makedirs(output_dir, exist_ok=True)

        kmc_cmd = f"kmc {merged_fastq} {os.path.join(sample_dir, current_fasta, current_fasta)} {kmc_dir} -t{thread} "

        if kmer_l:
            kmc_cmd += f" -k{kmer_l}"
        if min_kmer_coverage:
            kmc_cmd += f" -ci{min_kmer_coverage}"
        if max_kmer_coverage:
            kmc_cmd += f" -cs{max_kmer_coverage}"
        if kmer_times:
            kmc_cmd += f" -cx{kmer_times}"
        if memory:
            kmc_cmd += f" -m{memory}"

        subprocess.run(kmc_cmd, shell=True)

        histogram_path = os.path.join(sample_dir, f'{current_fasta}_sample.histo')
        kmc_tools_cmd = f"kmc_tools transform {os.path.join(sample_dir, f'{current_fasta}', f'{current_fasta}')} histogram {histogram_path} -cx{kmer_times}"
        
        subprocess.run(kmc_tools_cmd, shell=True)

        r_script = f"""
        kmer <- read.table('{histogram_path}', header=FALSE)
        if (nrow(kmer) > 0) {{
            kmer <- subset(kmer, V1 >= {min_freq} & V1 <= {max_freq})
            Frequency <- kmer$V1
            Number <- kmer$V2
            pdf('{os.path.join(sample_dir, 'kmer_plot.pdf')}')
            plot(Frequency, Number, type='l', col='blue', main='k-mer Frequency Plot', xlab='Frequency', ylab='Number of k-mers')
            dev.off()
        }} else {{
            cat("No data available for plotting.")
        }}
        """

        r_script_path = os.path.join(sample_dir, 'plot_kmer.R')
        with open(r_script_path, 'w') as f:
            f.write(r_script)

        subprocess.run(['Rscript', r_script_path], check=True)
        print(f"Plot saved to {os.path.join(sample_dir, 'kmer_plot.pdf')}")

        genomescope_cmd = f"Rscript {genomescope_path} -i {histogram_path} -k {kmer_l} -o {os.path.join(sample_dir, 'genome_feature')} -n {current_fasta} -m {kmer_times}"
        subprocess.run(genomescope_cmd, shell=True)




def custom_help():
    console = Console()

    intro = (
        "The 'Quality_control' script includes four functional modules: "
        "quality control of sequencing data, removal of host sequences, "
        "redundancy & distribution statistics, and k-mer frequency statistics & "
        "genome feature assessment. This script can be used for preliminary "
        "cleaning and evaluation of sequencing data."
    )

    examples_md = Markdown(
        "\n**Examples:**\n"
        "```\n"
        "Quality_control  -r /data  -out /results  -gi /hg38  -T_nnp kmer  -n 2048  -v 8  -kl 24  -m 66  -ci 1  -cs 256  -cx 10000 -min_freq 10  -max_freq 500\n"
        "Quality_control  -r /data  -out /results  -p 80\n"
        "```\n"
    )

    reads_md = Markdown(
        "**Input reads location (-r):** the directory contains paired-end FASTQ files (decompressed), e.g.\n"
        "```\n"
        "/xxxx/raw_data/A1/A1_1.fastq\n"
        "/xxxx/raw_data/A1/A1_2.fastq\n"
        "/xxxx/raw_data/A2/A2_1.fastq\n"
        "/xxxx/raw_data/A2/A2_2.fastq\n"
        "...\n"
        "```\n"
        "Therefore, `-r` should be set to `/xxxx/raw_data`.\n"
    )

    out_md = Markdown(
        "**Output directory (-out):** it will contain four subdirectories corresponding to four analyses."
    )

    console.print(Panel(intro, border_style="cyan", title="Quality_control", title_align="left"))
    console.print(examples_md)
    console.print(reads_md)
    console.print(out_md)
    console.print()

    # ---- Part 1: QC (fastp) --------------------------------------------------
    qc_tbl = Table(show_header=True, header_style="bold magenta")
    qc_tbl.add_column("Files", style="cyan", no_wrap=True)
    qc_tbl.add_column("Description", style="white")

    qc_tbl.add_row("A1_trim_R1.fastq.gz & A1_trim_R2.fastq.gz",
                   "Quality-controlled FASTQ files (adapter trimming, length filter, low-quality trimming).")
    qc_tbl.add_row("A1_trim-unpaired_R1.fastq.gz & A1_trim-unpaired_R2.fastq.gz",
                   "Filtered FASTQ files, kept for backup.")
    qc_tbl.add_row("A1_fastp_summary_before.tsv & A1_fastp_summary_after.tsv",
                   "Statistics on reads before and after quality control.")
    qc_tbl.add_row("A1_unpaired_summary.tsv", "Information on filtered reads.")
    qc_tbl.add_row("A1_fastp.json & A1_fastp.html", "Quality control result reports.")
    qc_tbl.add_row("~/origin_quality/A1_trim_R1_fastqc.html",
                   "Quality detection report of original reads before QC (FastQC HTML).")
    qc_tbl.add_row("~/filtered_quality/A1_trim-unpaired_R1_fastqc.html & ~/filtered_quality/A1_trim_R1_fastqc.html",
                   "Quality detection reports of reads after QC (FastQC HTML).")
    qc_tbl.add_row("~/origin_quality/*.zip & ~/filtered_quality/*.zip",
                   "After decompression, detailed quality detection result files. "
                   "Sub-directory 'Images' contains plots.")

    console.print(Panel(qc_tbl, border_style="magenta", title="1) QC (Basic Quality Control by Fastp)", title_align="left"))
    console.print()

    # ---- Part 2: rm_host -----------------------------------------------------
    rm_tbl = Table(show_header=True, header_style="bold green")
    rm_tbl.add_column("Files", style="cyan", no_wrap=True)
    rm_tbl.add_column("Description", style="white")

    rm_tbl.add_row("A1_host_removed_1.fq & A1_host_removed_2.fq",
                   "Paired-end FASTQ files after removing host reads.")
    rm_tbl.add_row("A1_unpaired_summary.tsv",
                   "Statistics on paired-end FASTQ files after host removal.")
    rm_tbl.add_row("(BAM/SAM intermediates)", "Intermediate results produced during the process.")
    console.print(Panel(rm_tbl, border_style="green", title="2) rm_host (Bowtie2, Samtools, Bedtools)", title_align="left"))
    console.print()

    # ---- Part 3: Nonpareil ---------------------------------------------------
    np_tbl = Table(show_header=True, header_style="bold yellow")
    np_tbl.add_column("Files", style="cyan", no_wrap=True)
    np_tbl.add_column("Description", style="white")

    np_tbl.add_row("A1_rmhost.fq", "Merged paired-end FASTQ file, higher quality for statistics.")
    np_tbl.add_row(".npc, .npl, .npo, .npa", "Nonpareil output files.")
    np_tbl.add_row("A1_nonpareil_index.tsv",
                   "Results from Nonpareil analysis (diversity, coverage, likelihood ratio, etc.).")
    np_tbl.add_row("A1_nonpareil_curves.pdf", "Visualization of Nonpareil results.")
    console.print(Panel(np_tbl, border_style="yellow", title="3) nonpareil (Redundancy & Distribution Statistics)", title_align="left"))
    console.print()

    # ---- Part 4: KMC ---------------------------------------------------------
    kmc_tbl = Table(show_header=True, header_style="bold blue")
    kmc_tbl.add_column("Files", style="cyan", no_wrap=True)
    kmc_tbl.add_column("Description", style="white")

    kmc_tbl.add_row("A1_sample.histo & kmer_plot.pdf",
                    "K-mer frequency distribution table and visualization.")
    kmc_tbl.add_row("/genome_feature/*", "Genome feature assessment and visualization.")
    console.print(Panel(kmc_tbl, border_style="blue", title="4) KMC (K-mer Frequency & Genome Feature Assessment)", title_align="left"))
    console.print(Panel("All four sections can be analyzed using preset parameters, and users can adjust parameters according to their data needs.",
                        border_style="cyan"))
    console.print(Rule())

    # ====================== Detailed parameters ===============================
    console.print("[bold]Detailed parameters[/bold]\n")

    # --- Usage & Global
    console.print(Panel("Usage: Quality_control [OPTIONS] INPUT OUTPUT DATABASE",
                        border_style="cyan", title="Usage & Global parameters", title_align="left"))

    g_tbl = Table(show_header=False, box=None, pad_edge=False)
    g_tbl.add_column("Flag", style="bold cyan", no_wrap=True)
    g_tbl.add_column("Description", style="white")
    g_tbl.add_row("-cf, --config", "DPath to config.yml for env/tools/dbs.")
    g_tbl.add_row("-r, --raw_data", "Directory containing subdirectories for single-sample data.")
    g_tbl.add_row("-out, --out_dir", "Output directory for quality control and assembly.")
    g_tbl.add_row("-p, --thread", "Number of parallel threads (default: 80).")
    console.print(g_tbl)
    console.print()

    # --- Part1: fastp options
    console.print(Panel("Part1: Basic quality control by Fastp", border_style="magenta"))

    t1 = Table(show_header=False, box=None, pad_edge=False)
    t1.add_column("Flag", style="bold cyan", no_wrap=True)
    t1.add_column("Description", style="white")

    # Quality Filtering
    console.print("[bold]Quality Filtering:[/bold]")
    t1.add_row("-q, --qualified_quality", "Set the base quality value (default: 15).")
    t1.add_row("-u, --unqualified_percent", "Proportion of unqualified bases allowed (default: 40).")
    console.print(t1)

    # Length Filtering
    t1a = Table(show_header=False, box=None, pad_edge=False)
    t1a.add_column("Flag", style="bold cyan", no_wrap=True)
    t1a.add_column("Description", style="white")
    console.print()
    console.print("[bold]Length Filtering:[/bold]")
    t1a.add_row("-l_min, --length_required", "Minimum length of reads (default: 15).")
    t1a.add_row("-l_max, --length_limit", "Maximum length of reads (default: no limit).")
    console.print(t1a)

    # Low Complexity
    t1b = Table(show_header=False, box=None, pad_edge=False)
    t1b.add_column("Flag", style="bold cyan", no_wrap=True)
    t1b.add_column("Description", style="white")
    console.print()
    console.print("[bold]Low Complexity Filtering:[/bold]")
    t1b.add_row("-Y, --complexity_threshold", "Complexity-filtering threshold for reads (default: 30).")
    console.print(t1b)

    # Adapter
    t1c = Table(show_header=False, box=None, pad_edge=False)
    t1c.add_column("Flag", style="bold cyan", no_wrap=True)
    t1c.add_column("Description", style="white")
    console.print()
    console.print("[bold]Adapter Trimming:[/bold]")
    t1c.add_row("-a, --adapter", "Path of adapter fasta file.")
    console.print(t1c)

    # PE Correction
    t1d = Table(show_header=False, box=None, pad_edge=False)
    t1d.add_column("Flag", style="bold cyan", no_wrap=True)
    t1d.add_column("Description", style="white")
    console.print()
    console.print("[bold]Paired-End Data Correction:[/bold]")
    t1d.add_row("-c, --correction", "Perform clipping and correction on paired-end (PE) data.")
    console.print(t1d)

    # R1 trimming
    t1e = Table(show_header=False, box=None, pad_edge=False)
    t1e.add_column("Flag", style="bold cyan", no_wrap=True)
    t1e.add_column("Description", style="white")
    console.print()
    console.print("[bold]Trimming Parameters for R1:[/bold]")
    t1e.add_row("-f, --trim_front1", "Remove base pairs from the start of R1 (default: 0).")
    t1e.add_row("-t, --trim_tail1", "Remove base pairs from the end of R1 (default: 0).")
    t1e.add_row("-b, --max_len1", "Maximum length threshold for R1 (default: no limit).")
    console.print(t1e)

    # R2 trimming
    t1f = Table(show_header=False, box=None, pad_edge=False)
    t1f.add_column("Flag", style="bold cyan", no_wrap=True)
    t1f.add_column("Description", style="white")
    console.print()
    console.print("[bold]Trimming Parameters for R2:[/bold]")
    t1f.add_row("-F, --trim_front2", "Remove base pairs from the start of R2 (default: same as R1).")
    t1f.add_row("-T, --trim_tail2", "Remove base pairs from the end of R2 (default: same as R1).")
    t1f.add_row("-B, --max_len2", "Maximum length for R2 (default: no limit).")

    console.print(t1f)
    console.print()

    # --- Part2: host removal
    console.print(Panel("Part2: Remove host genome by Bowtie2, Samtools and bedtools", border_style="green"))
    t2 = Table(show_header=False, box=None, pad_edge=False)
    t2.add_column("Flag", style="bold cyan", no_wrap=True)
    t2.add_column("Description", style="white")
    t2.add_row("-gf, --genome_fasta", "Path to the host genome FASTA file.")
    t2.add_row("-gi, --genome_index", "Index file for the host genome.")
    console.print(t2)
    console.print()

    # --- Part3: Nonpareil
    console.print(Panel("Part3: Redundancy & Distribution Statistics", border_style="yellow"))
    t3 = Table(show_header=False, box=None, pad_edge=False)
    t3.add_column("Flag", style="bold cyan", no_wrap=True)
    t3.add_column("Description", style="white")
    console.print("[bold]Mandatory parameters:[/bold]")
    t3.add_row("-min_cut, --min_reads_cut", "Remove short reads before redundancy and distribution evaluation. (default: 15)")
    t3.add_row("-T_nnp, --nonpareil_mth", "Method for redundancy analysis ('kmer' or 'alignment').")
    console.print(t3)

    t3c = Table(show_header=False, box=None, pad_edge=False)
    t3c.add_column("Flag", style="bold cyan", no_wrap=True)
    t3c.add_column("Description", style="white")
    console.print()
    console.print("[bold]Common parameters:[/bold]")
    t3c.add_row("-X, --max_reads", "Maximum reads for analysis (default: 10000 for 'kmer', 1000 for 'alignment').")
    t3c.add_row("-k, --kmer_length", "Length of kmer (default: 24).")
    t3c.add_row("-n, --num_subsamples", "Number of subsamples per point (default: 1024).")
    t3c.add_row("-L, --min_overlap", "Minimum overlap percentage for alignment regions (default: 50.0).")
    t3c.add_row("-R, --max_ram", "Maximum RAM usage in MiB (default: 1024).")
    t3c.add_row("-v, --verbose", "Verbosity level (default: 7).")
    t3c.add_row("-r_seed, --random_seed", "Random seed for reproducibility (only for 'alignment').")
    console.print(t3c)
    console.print()

    # --- Part4: KMC
    console.print(Panel("Part4: K-mer frequency statistics & Genome feature assessment", border_style="blue"))
    t4 = Table(show_header=False, box=None, pad_edge=False)
    t4.add_column("Flag", style="bold cyan", no_wrap=True)
    t4.add_column("Description", style="white")
    t4.add_row("-kl, --kmer_l", "Specify k-mer length (range: 1 to 256, default: 25).")
    t4.add_row("-m, --memory", "Set maximum memory usage (GB, range: 1 to 1024, default: 12).")
    t4.add_row("-ci, --min_kmer_coverage", "Exclude k-mers with coverage below this value (default: 2).")
    t4.add_row("-cs, --max_kmer_coverage", "Set maximum value for k-mer counters (default: 255).")
    t4.add_row("-cx, --kmer_times", "Exclude k-mers that appear more than this value (default: 10000).")
    t4.add_row("-min_freq, --min_freq", "Minimum k-mer length when plot k-mer frequency (default: 5).")
    t4.add_row("-max_freq, --max_freq", "Maximum k-mer length when plot k-mer frequency (default: 500).")
    console.print(t4)

    console.print()
    console.print(Rule(style="dim"))



class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="Upstream of metagenomics (quality control and contigs assembly).")
    parser.add_argument('-cf','--config', type=str, help='Path to config.yml for env/tools/dbs')
    parser.add_argument("-r", "--raw_data", type=str, required=True)
    parser.add_argument("-out", "--out_dir", type=str, required=True)
    parser.add_argument("-p", "--thread", type=int, default=80)
    parser.add_argument("-q", "--qualified_quality", type=int, default=15)
    parser.add_argument("-u", "--unqualified_percent", type=int, default=40)
    parser.add_argument("-l_min", "--length_required", type=int, default=15)
    parser.add_argument("-l_max", "--length_limit", type=int, default=0)
    parser.add_argument("-Y", "--complexity_threshold", type=int, default=30)
    parser.add_argument("-a", "--adapter", type=str)
    parser.add_argument("-c", "--correction", action='store_true')
    parser.add_argument("-f", "--trim_front1", type=int, default=0)
    parser.add_argument("-t", "--trim_tail1", type=int, default=0)
    parser.add_argument("-b", "--max_len1", type=int, default=0)
    parser.add_argument("-F", "--trim_front2", type=int, default=0)
    parser.add_argument("-T", "--trim_tail2", type=int, default=0)
    parser.add_argument("-B", "--max_len2", type=int, default=0)
    parser.add_argument("-gf", "--genome_fasta", type=str)
    parser.add_argument("-gi", "--genome_index", type=str)
    parser.add_argument("-min_cut", "--min_reads_cut", type=int, default=15)
    parser.add_argument("-T_nnp", "--nonpareil_mth", type=str, choices=['kmer', 'alignment'], default='kmer')
    parser.add_argument("-X", "--max_reads", type=int, default=10000)
    parser.add_argument("-k", "--kmer_length", type=int, default=24)   
    parser.add_argument("-n", "--num_subsamples", type=int, default=1024)
    parser.add_argument("-L", "--min_overlap", type=float, default=50.0)
    parser.add_argument("-R", "--max_ram", type=int, default=1024)
    parser.add_argument("-v", "--verbose", type=int, default=7)
    parser.add_argument("-r_seed", "--random_seed", type=int)
    parser.add_argument("-kl", "--kmer_l", type=int, default=25)
    parser.add_argument("-m", "--memory", type=int, default=12)
    parser.add_argument("-ci", "--min_kmer_coverage", type=int, default=2)
    parser.add_argument("-cs", "--max_kmer_coverage", type=int, default=255)
    parser.add_argument("-cx", "--kmer_times", type=int, default=10000)
    parser.add_argument("-min_freq", "--min_freq", type=int, default=5)
    parser.add_argument("-max_freq", "--max_freq", type=int, default=500)
    args = parser.parse_args()

    cfg = {}
    if args.config:
        with open(args.config, 'r') as f:
            raw = yaml.safe_load(f) or {}
            cfg.update(raw)
    global CONFIG
    CONFIG = cfg
    vmp_env = cfg.get('VMP_env')
    if vmp_env:
        vmp_env = str(vmp_env).strip().strip("'\"")
        env_bin = os.path.join(vmp_env, 'bin')
        if os.path.isdir(env_bin):
            os.environ['PATH'] = env_bin + os.pathsep + os.environ.get('PATH','')
        else:
            print(f"[warn] VMP_env bin not found: {env_bin} (continuing)")
    if (getattr(args, 'genome_fasta', None) is None or getattr(args, 'genome_index', None) is None) and cfg.get('host_genome'):
        host_dir = str(cfg.get('host_genome')).strip().strip("'\"")
        default_fa = os.path.join(host_dir, 'hg38.fa')
        # put bowtie2 index under output dir to avoid permission issues
        default_index_dir = os.path.join(args.out_dir, 'rm_host', 'host_index')
        os.makedirs(default_index_dir, exist_ok=True)
        default_index = os.path.join(default_index_dir, 'hg38')
        if getattr(args, 'genome_fasta', None) is None:
            args.genome_fasta = default_fa
        if getattr(args, 'genome_index', None) is None:
            args.genome_index = default_index

    cfg = {}
    if args.config:
        with open(args.config, 'r') as f:
            raw = yaml.safe_load(f) or {}
            cfg.update(raw)
    vmp_env = cfg.get('VMP_env')
    if vmp_env:
        vmp_env = str(vmp_env).strip().strip("'\"")
        env_bin = os.path.join(vmp_env, 'bin')
        if os.path.isdir(env_bin):
            os.environ['PATH'] = env_bin + os.pathsep + os.environ.get('PATH','')
        else:
            print(f"[warn] VMP_env bin not found: {env_bin} (continuing)")

    if (getattr(args, 'genome_fasta', None) is None or getattr(args, 'genome_index', None) is None) and cfg.get('host_genome'):
        host_dir = str(cfg.get('host_genome')).strip().strip("'\"")
        default_fa = os.path.join(host_dir, 'hg38.fa')
        default_index_dir = os.path.join(args.out_dir, 'rm_host', 'host_index')
        os.makedirs(default_index_dir, exist_ok=True)
        default_index = os.path.join(default_index_dir, 'hg38')
        if getattr(args, 'genome_fasta', None) is None:
            args.genome_fasta = default_fa
        if getattr(args, 'genome_index', None) is None:
            args.genome_index = default_index


    fastq_list = Getfastq_list(args.raw_data)
    fastp(args.out_dir, args.raw_data, fastq_list, args.qualified_quality, args.unqualified_percent, args.length_required, args.length_limit, args.complexity_threshold, args.adapter, args.correction, args.trim_front1, args.trim_tail1, args.max_len1, args.trim_front2, args.trim_tail2, args.max_len2, args.thread)
    rm_host(args.out_dir, fastq_list, args.genome_fasta, args.genome_index, args.thread, args.raw_data)
    Nonpareil(args.out_dir, fastq_list, args.thread, args.min_reads_cut, args.nonpareil_mth, args.max_reads, args.kmer_length, args.num_subsamples, args.min_overlap, args.max_ram, args.verbose, args.random_seed)
    KMC(args.out_dir, fastq_list, args.thread, args.kmer_l, args.min_kmer_coverage, args.max_kmer_coverage, args.kmer_times, args.memory, args.min_freq, args.max_freq)


if __name__ == "__main__":
    main()

