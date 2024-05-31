import os
import subprocess
import wx
import threading
from Bio import SeqIO
from Bio.Blast import NCBIXML
import io

def run_blast(selected_query_path, selected_cry_genes_path, evalue_threshold, word_size, gap_open, gap_extend, result_handler):
    try:
        query_sequences = list(SeqIO.parse(selected_query_path, "fasta"))

        if not query_sequences:
            wx.MessageBox("No sequences found in the selected query file.", "Error", wx.OK | wx.ICON_ERROR)
            return

        result_handler.result_text.SetValue("")  # Clear previous results
        result_handler.status_label.SetLabel("Status: Running BLASTn...")
        result_handler.run_blast_button.Disable()  # Disable button during BLAST

        def process_blast_results(query_sequences, cry_genes_path, evalue_threshold, word_size, gap_open, gap_extend, result_handler):
            results = ""
            evalues = []

            # Gather all .fasta files in the cry_genes_path directory
            fasta_files = [os.path.join(cry_genes_path, f) for f in os.listdir(cry_genes_path) if f.endswith(".fasta")]

            # Combine all .fasta files into a single database file with modified titles
            combined_fasta_path = os.path.join(cry_genes_path, "combined_cry_genes.fasta")
            with open(combined_fasta_path, 'w') as combined_fasta:
                for fasta_file in fasta_files:
                    gene_name = os.path.basename(fasta_file).replace('.fasta', '')
                    with open(fasta_file, 'r') as f:
                        for line in f:
                            if line.startswith(">"):
                                combined_fasta.write(f">{gene_name}\n")
                            else:
                                combined_fasta.write(line)

            # Create the BLAST database
            db_name = os.path.join(cry_genes_path, "cry_gene_db")
            makeblastdb_cmd = ['makeblastdb', '-in', combined_fasta_path, '-dbtype', 'nucl', '-out', db_name]
            subprocess.run(makeblastdb_cmd, check=True)

            for seq_record in query_sequences:
                query_gene = str(seq_record.seq)
                query_gene_name = os.path.basename(selected_query_path).replace('.fasta', '')
                query_file = os.path.join(cry_genes_path, "query.fasta")
                with open(query_file, 'w') as f:
                    f.write(f">{query_gene_name}\n{query_gene}\n")

                blastn_cmd = [
                    'blastn',
                    '-query', query_file,
                    '-db', db_name,
                    '-evalue', str(evalue_threshold),
                    '-word_size', str(word_size),
                    '-gapopen', str(gap_open),
                    '-gapextend', str(gap_extend),
                    '-outfmt', '5'
                ]
                result_handle = subprocess.run(blastn_cmd, capture_output=True, text=True, check=True)
                blast_output = result_handle.stdout

                # Use StringIO to convert the string output to a file-like object
                blast_output_io = io.StringIO(blast_output)
                blast_records = NCBIXML.parse(blast_output_io)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        cry_gene_name = alignment.hit_def.split()[0]  # Use the filename part from the modified title
                        for hsp in alignment.hsps:
                            evalue = hsp.expect
                            evalues.append(evalue)
                            results += f"Query Sequence: {query_gene_name}\n"
                            results += f"Cry Gene: {cry_gene_name}\n"
                            results += f"Length: {alignment.length}\n"
                            results += f"E-value: {evalue}\n"
                            results += f"Score: {hsp.score}\n\n"

                            # Check if Cry gene is detected based on score and e-value criteria
                            if hsp.score > 400 and evalue < 1e-10:
                                results += f"Plant is likely a GMO (contains {cry_gene_name} matched with a score of {hsp.score} and e-value of {evalue})\n"

            result_handler.display_results(results, evalues)

        blast_thread = threading.Thread(target=process_blast_results,
                                        args=(query_sequences, selected_cry_genes_path, evalue_threshold, word_size, gap_open, gap_extend, result_handler))
        blast_thread.start()

    except Exception as e:
        result_handler.handle_error(str(e))
