import os
import subprocess
import pandas as pd
import json
import argparse
import shlex
import shutil
# from termcolor import colored

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
        # print(f"Sample {current_fasta} has completed quality control.")

        stat_unpaired_reads(os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R1.fastq.gz'), current_fasta, sample_dir)
        stat_unpaired_reads(os.path.join(sample_dir, f'{current_fasta}_trim-unpaired_R2.fastq.gz'), current_fasta, sample_dir)

        generate_fastp_summary(os.path.join(sample_dir, f'{current_fasta}_fastp.json'), current_fasta, sample_dir)

        # # FastQC
        # origin_quality_dir = os.path.join(sample_dir, "origin_quality")
        # filtered_quality_dir = os.path.join(sample_dir, "filtered_quality")
        # os.makedirs(origin_quality_dir, exist_ok=True)
        # os.makedirs(filtered_quality_dir, exist_ok=True)
        # ##### quality evaluation
        # fastqc_filtered_command = f"fastqc -t {thread} {os.path.join(sample_dir, '*.fastq.gz')} -o {filtered_quality_dir}"
        # fastqc_orgin_command = f"fastqc -t {thread} {os.path.join(raw_data, current_fasta, f'*.fastq')} -o {origin_quality_dir}"
        # os.system(fastqc_filtered_command)
        # os.system(fastqc_orgin_command)

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

        # Convert SAM to BAM    
        # samtools view -@ 80 -bSh /opt/user/basspoom/QC_new/Up-stream/rm_host/A1/A1_host_mapped.sam > /opt/user/basspoom/QC_new/Up-stream/rm_host/A1/A1_host_mapped.bam
        # samtools view -@ 80 -bSh "/opt/user/basspoom/QC_new/Up-stream/rm_host/A1/A1_host_mapped.sam" > "/opt/user/basspoom/QC_new/Up-stream/rm_host/A1/A1_host_mapped.bam"

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
    YELLOW = "\033[93m"  # Yellow
    GREEN = "\033[92m"   # Green
    BLUE = "\033[94m"    # Blue
    PURPLE = "\033[95m"  # Purple
    RED = "\033[91m"     # Red
    RESET = "\033[0m"    # Reset to default color

    print("\n" + RED + "The 'Quality_control' script includes four functional modules: quality control of sequencing data, removal of host sequences, redundancy & distribution statistics, and k-mer frequency statistics & genome feature assessment. This script can be used for preliminary cleaning and evaluation of sequencing data." + RESET)

    print("\n" + PURPLE + "Examples:" + RESET)
    print("  Quality_control  -r /data  -out /results  -gi /hg38  -T_nnp kmer  -n 2048  -v 8  -kl 24  -m 66  -ci 1  -cs 256  -cx 10000 -min_freq 10  -max_freq 500")
    print("  Quality_control  -r /data  -out /results  -p 80")

    print("\n" + GREEN + "The -r parameter points to the directory containing the paired-end sequencing FASTQ files (please decompress first), for example:" + RESET)
    print("  " + BLUE + "/xxxx/raw_data/A1/A1_1.fastq" + RESET)
    print("  " + BLUE + "/xxxx/raw_data/A1/A1_2.fastq" + RESET)
    print("  " + BLUE + "/xxxx/raw_data/A2/A2_1.fastq" + RESET)
    print("  " + BLUE + "/xxxx/raw_data/A2/A2_2.fastq" + RESET)
    print("  ......")
    print("Therefore, -r should be set to /xxxx/raw_data.")

    print("\n" + GREEN + "The -out parameter specifies the output directory, which will contain four subdirectories corresponding to four analyses:" + RESET)

    print("\n" + GREEN + "1. QC (Part 1: Basic Quality Control by Fastp)" + RESET)
    print("  " + BLUE + "A1_trim_R1.fastq.gz" + RESET + " & " + BLUE + "A1_trim_R2.fastq.gz" + RESET + ": Quality-controlled FASTQ files, processed to remove adapters, filter by length, and trim low-quality reads.")
    print("  " + BLUE + "A1_trim-unpaired_R1.fastq.gz" + RESET + " & " + BLUE + "A1_trim-unpaired_R2.fastq.gz" + RESET + ": Filtered FASTQ files, kept for backup.")
    print("  " + BLUE + "A1_fastp_summary_before.tsv" + RESET + " & " + BLUE + "A1_fastp_summary_after.tsv" + RESET + ": Statistics on reads before and after quality control.")
    print("  " + BLUE + "A1_unpaired_summary.tsv" + RESET + ": Information on filtered reads.")
    print("  " + BLUE + "A1_fastp.json" + RESET + " & " + BLUE + "A1_fastp.html" + RESET + ": Quality control result reports.")
    # print("\n" + GREEN + "Quality Reports" + RESET)
    print("  " + BLUE + "~/origin_quality/A1_trim_R1_fastqc.html" + RESET + ": Quality detection report of original reads before quality control, in HTML format, generated by FastQC.")
    print("  " + BLUE + "~/filtered_quality/A1_trim-unpaired_R1_fastqc.html" + RESET + " & " + BLUE + "~/filtered_quality/A1_trim_R1_fastqc.html" + RESET + ": Quality detection reports of reads after quality control, in HTML format, generated by FastQC.")
    print("  " + BLUE + "~/origin_quality/*.zip" + RESET + " & " + BLUE + "~/filtered_quality/*.zip" + RESET + ": After decompression, detailed quality detection result files are obtained. The sub - directory 'Images' contains image displays of quality results.")
    
    print("\n" + GREEN + "2. rm_host (Part 2: Remove Host Genome by Bowtie2, Samtools, and Bedtools)" + RESET)
    print("  " + BLUE + "A1_host_removed_1.fq" + RESET + " & " + BLUE + "A1_host_removed_2.fq" + RESET + ": Paired-end FASTQ files after removing host reads.")
    print("  " + BLUE + "A1_unpaired_summary.tsv" + RESET + ": Statistics on paired-end FASTQ files after host removal.")
    print("  (Other BAM and SAM files are intermediate results produced during the process.)")

    print("\n" + GREEN + "3. nonpareil (Part 3: Redundancy & Distribution Statistics)" + RESET)
    print("  " + BLUE + "A1_rmhost.fq" + RESET + ": Merged paired-end FASTQ file, higher quality for statistics.")
    print("  " + BLUE + ".npc" + RESET + ", " + BLUE + ".npl" + RESET + ", " + BLUE + ".npo" + RESET + ", " + BLUE + ".npa" + RESET + ": Nonpareil output files.")
    print("  " + BLUE + "A1_nonpareil_index.tsv" + RESET + ": Results from Nonpareil analysis, including diversity, coverage, likelihood ratio, etc.")
    print("  " + BLUE + "A1_nonpareil_curves.pdf" + RESET + ": Visualization of Nonpareil results.")

    print("\n" + GREEN + "4. KMC (Part 4: K-mer Frequency Statistics & Genome Feature Assessment)" + RESET)
    print("  " + BLUE + "A1_sample.histo" + RESET + " & " + BLUE + "kmer_plot.pdf" + RESET + ": K-mer frequency distribution table and visualization.")
    print("  " + BLUE + "/genome_feature/*" + RESET + ": Genome feature assessment and visualization.")

    print("\n" + GREEN + "All four sections can be analyzed using preset parameters, and users can adjust parameters according to their data needs." + RESET)


    print("\n" + PURPLE + "Detailed parameters")
    # Print usage and global parameters
    print(" " + GREEN + "=" * 50 + " " + RESET + YELLOW + "Usage & Global parameters" + RESET + " " + GREEN + "=" * 50 + RESET)#########################

    print("\n" + YELLOW + "Usage: Quality_control [OPTIONS] INPUT OUTPUT DATABASE" + RESET + "\n")

    print(YELLOW + "Global parameters:" + RESET)
    print(f"  {'-r, --raw_data':<40} Directory containing subdirectories for single-sample data.")
    print(f"  {'-out, --out_dir':<40} Output directory for quality control and assembly.")
    # print(f"  {'-qc_env, --quality_control_env':<40} QC's environment location (e.g., '/opt/user/xxx/tools/anaconda3/envs/QC').")
    print(f"  {'-p, --thread':<40} Number of parallel threads (default: 80).")

    print("\n" + GREEN + "=" * 44 + " " + RESET + YELLOW + "Part1: Basic quality control by Fastp" + RESET + " " + GREEN + "=" * 44 + RESET)################################

    print("\n" + YELLOW + "Quality Filtering:" + RESET)
    print(f"  {'-q, --qualified_quality':<40} Set the base quality value (default: 15).")
    print(f"  {'-u, --unqualified_percent':<40} Proportion of unqualified bases allowed (default: 40).")
    
    print("\n" + YELLOW + "Length Filtering:" + RESET)
    print(f"  {'-l_min, --length_required':<40} Minimum length of reads (default: 15).")
    print(f"  {'-l_max, --length_limit':<40} Maximum length of reads (default: no limit).")
    
    print("\n" + YELLOW + "Low Complexity Filtering:" + RESET)
    print(f"  {'-Y, --complexity_threshold':<40} Complexity-filtering threshold for reads (default: 30).")
    
    print("\n" + YELLOW + "Adapter Trimming:" + RESET)
    print(f"  {'-a, --adapter':<40} Path of adapter fasta file.")
    
    print("\n" + YELLOW + "Paired-End Data Correction:" + RESET)
    print(f"  {'-c, --correction':<40} Perform clipping and correction on paired-end (PE) data.")
    
    print("\n" + YELLOW + "Trimming Parameters for R1:" + RESET)
    print(f"  {'-f, --trim_front1':<40} Remove base pairs from the start of R1 (default: 0).")
    print(f"  {'-t, --trim_tail1':<40} Remove base pairs from the end of R1 (default: 0).")
    print(f"  {'-b, --max_len1':<40} Maximum length threshold for R1 (default: no limit).")
    
    print("\n" + YELLOW + "Trimming Parameters for R2:" + RESET)
    print(f"  {'-F, --trim_front2':<40} Remove base pairs from the start of R2 (default: same as R1).")
    print(f"  {'-T, --trim_tail2':<40} Remove base pairs from the end of R2 (default: same as R1).")
    print(f"  {'-B, --max_len2':<40} Maximum length for R2 (default: no limit).")

    print("\n" + GREEN + "=" * 33 + " " + RESET + YELLOW + "Part2: Remove host genome by Bowtie2, Samtools and bedtools" + RESET + " " + GREEN + "=" * 33 + RESET)###############################

    print("\n" + YELLOW + "Host Genome Parameters:" + RESET)
    print(f"  {'-gf, --genome_fasta':<40} Path to the host genome FASTA file.")
    print(f"  {'-gi, --genome_index':<40} Index file for the host genome.")
    
    print("\n" + GREEN + "=" * 41 + " " + RESET + YELLOW + "Part3: Redundancy & Distribution Statistics" + RESET + " " + GREEN + "=" * 41 + RESET)####################################

    print("\n" + YELLOW + "Mandatory parameters:" + RESET)
    print(f"  {'-min_cut, --min_reads_cut':<40} Remove short reads before redundancy and distribution evaluation. (default: 15)")  
    print(f"  {'-T_nnp, --nonpareil_mth':<40} Method for redundancy analysis ('kmer' or 'alignment').")

    print("\n" + YELLOW + "Common parameters:" + RESET)
    print(f"  {'-X, --max_reads':<40} Maximum reads for analysis (default: 10000 for 'kmer', 1000 for 'alignment').")
    print(f"  {'-k, --kmer_length':<40} Length of kmer (default: 24).")
    print(f"  {'-n, --num_subsamples':<40} Number of subsamples per point (default: 1024).")
    print(f"  {'-L, --min_overlap':<40} Minimum overlap percentage for alignment regions (default: 50.0).")
    print(f"  {'-R, --max_ram':<40} Maximum RAM usage in MiB (default: 1024).")
    print(f"  {'-v, --verbose':<40} Verbosity level (default: 7).")
    print(f"  {'-r_seed, --random_seed':<40} Random seed for reproducibility (only for 'alignment').")

    print("\n" + GREEN + "=" * 32 + " " + RESET + YELLOW + "Part4: K-mer frequency statistics & Genome feature assessment" + RESET + " " + GREEN + "=" * 32 + RESET)##############################

    print("\n" + YELLOW + "Frequency statistics of k-mer & Genome feature assessment:" + RESET)
    print(f"  {'-kl, --kmer_l':<40} Specify k-mer length (range: 1 to 256, default: 25).")
    print(f"  {'-m, --memory':<40} Set maximum memory usage (GB, range: 1 to 1024, default: 12).")
    print(f"  {'-ci, --min_kmer_coverage':<40} Exclude k-mers with coverage below this value (default: 2).")
    print(f"  {'-cs, --max_kmer_coverage':<40} Set maximum value for k-mer counters (default: 255).")
    print(f"  {'-cx, --kmer_times':<40} Exclude k-mers that appear more than this value (default: 10000).")
    print(f"  {'-min_freq, --min_freq':<40} Minimum k-mer length when plot k-mer frequency (default: 5).")
    print(f"  {'-max_freq, --max_freq':<40} Maximum k-mer length when plot k-mer frequency (default: 500).")



class CustomArgumentParser(argparse.ArgumentParser):
    def print_help(self, file=None):
        custom_help()
        self.exit()

def main():
    parser = CustomArgumentParser(description="Upstream of metagenomics (quality control and contigs assembly).")
    parser.add_argument("-r", "--raw_data", type=str, required=True, help="Directory containing subdirectories for single-sample data.")
    parser.add_argument("-out", "--out_dir", type=str, required=True, help="Output directory for quality control and assembly.")
    parser.add_argument("-p", "--thread", type=int, default=80, help="Number of parallel threads (default: 80)")

    parser.add_argument("-q", "--qualified_quality", type=int, default=15, help="Base quality value (default: 15).")
    parser.add_argument("-u", "--unqualified_percent", type=int, default=40, help="Proportion of unqualified bases (default: 40).")
    parser.add_argument("-l_min", "--length_required", type=int, default=15, help="Minimum length of reads (default: 15).")
    parser.add_argument("-l_max", "--length_limit", type=int, default=0, help="Maximum length of reads (default: no limit).")
    parser.add_argument("-Y", "--complexity_threshold", type=int, default=30, help="Complexity-filtering threshold (default: 30).")
    parser.add_argument("-a", "--adapter", type=str, help="Path of adapter fasta file.")
    parser.add_argument("-c", "--correction", action='store_true', help="Perform clipping and correction on paired-end (PE) data.")
    parser.add_argument("-f", "--trim_front1", type=int, default=0, help="Remove base pairs from start of R1 (default: 0).")
    parser.add_argument("-t", "--trim_tail1", type=int, default=0, help="Remove base pairs from end of R1 (default: 0).")
    parser.add_argument("-b", "--max_len1", type=int, default=0, help="Maximum length for R1 (default: 0 == no limit).")
    parser.add_argument("-F", "--trim_front2", type=int, default=0, help="Remove base pairs from start of R2 (default: same as R1).")
    parser.add_argument("-T", "--trim_tail2", type=int, default=0, help="Remove base pairs from end of R2 (default: same as R1).")
    parser.add_argument("-B", "--max_len2", type=int, default=0, help="Maximum length for R2 (default: 0 == no limit).")

    # 去宿主部分参数
    parser.add_argument("-gf", "--genome_fasta", type=str, help="Path of reference genome.fasta.")
    parser.add_argument("-gi", "--genome_index", type=str, help="Directory containing bowtie2 index for reference genome.")

    # 冗余度&分布统计
    parser.add_argument("-min_cut", "--min_reads_cut", type=int, default=15, help='Remove short reads before redundancy and distribution evaluation. (default: 15)')
    parser.add_argument("-T_nnp", "--nonpareil_mth", type=str, choices=['kmer', 'alignment'], default='kmer', help='"kmer" or "alignment" (default == kmer)')
    parser.add_argument("-X", "--max_reads", type=int, default=10000, help='Default "kmer"=10000 or "alignment"=1000')
    parser.add_argument("-k", "--kmer_length", type=int, default=24, help='Length of kmer, default 24.')   
    parser.add_argument("-n", "--num_subsamples", type=int, default=1024, help='Number of subsamples per point, default 1024.')
    parser.add_argument("-L", "--min_overlap", type=float, default=50.0, help='Minimum overlap percentage for alignment regions, default 50.')
    parser.add_argument("-R", "--max_ram", type=int, default=1024, help='Maximum RAM usage (MiB), default 1024.')
    parser.add_argument("-v", "--verbose", type=int, default=7, help='Verbosity level, default 7.')
    parser.add_argument("-r_seed", "--random_seed", type=int, help='Random seed for reproducibility, only implemented for "alignment" method.')

    #  k-mer频数统计&基因组特征评估  冗余度&分布统计 
    parser.add_argument("-kl", "--kmer_l", type=int, default=25, help="Specify k-mer length (range: 1 to 256, default: 25).")
    parser.add_argument("-m", "--memory", type=int, default=12, help="Set maximum memory usage (GB, range: 1 to 1024, default: 12).")
    parser.add_argument("-ci", "--min_kmer_coverage", type=int, default=2, help="Exclude k-mers with coverage below this value (default: 2).")
    parser.add_argument("-cs", "--max_kmer_coverage", type=int, default=255, help="Set maximum value for k-mer counters (default: 255).")
    parser.add_argument("-cx", "--kmer_times", type=int, default=10000, help="Exclude k-mers that appear more than this value (default: 10000).")
    parser.add_argument("-min_freq", "--min_freq", type=int, default=5, help="Minimum k-mer length when plot k-mer Frequency (default: 5).")
    parser.add_argument("-max_freq", "--max_freq", type=int, default=500, help="Maxmum k-mer length when plot k-mer Frequency (default: 500).")

    args = parser.parse_args()

    fastq_list = Getfastq_list(args.raw_data)
    fastp(args.out_dir, args.raw_data, fastq_list, args.qualified_quality, args.unqualified_percent, args.length_required, args.length_limit, args.complexity_threshold, args.adapter, args.correction, args.trim_front1, args.trim_tail1, args.max_len1, args.trim_front2, args.trim_tail2, args.max_len2, args.thread)
    rm_host(args.out_dir, fastq_list, args.genome_fasta, args.genome_index, args.thread, args.raw_data)
    Nonpareil(args.out_dir, fastq_list, args.thread, args.min_reads_cut, args.nonpareil_mth, args.max_reads, args.kmer_length, args.num_subsamples, args.min_overlap, args.max_ram, args.verbose, args.random_seed)
    KMC(args.out_dir, fastq_list, args.thread, args.kmer_l, args.min_kmer_coverage, args.max_kmer_coverage, args.kmer_times, args.memory, args.min_freq, args.max_freq)


if __name__ == "__main__":
    main()
