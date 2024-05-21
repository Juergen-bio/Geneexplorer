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

            # Combine all .fasta files into a single database file
            combined_fasta_path = os.path.join(cry_genes_path, "combined_cry_genes.fasta")
            with open(combined_fasta_path, 'w') as combined_fasta:
                for fasta_file in fasta_files:
                    with open(fasta_file, 'r') as f:
                        combined_fasta.write(f.read())

            # Create the BLAST database
            db_name = os.path.join(cry_genes_path, "cry_gene_db")
            makeblastdb_cmd = ['makeblastdb', '-in', combined_fasta_path, '-dbtype', 'nucl', '-out', db_name]
            subprocess.run(makeblastdb_cmd, check=True)

            for seq_record in query_sequences:
                query_gene = str(seq_record.seq)
                query_file = os.path.join(cry_genes_path, "query.fasta")
                with open(query_file, 'w') as f:
                    f.write(f">{seq_record.id}\n{query_gene}\n")

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

                cry_gene_found = False

                # Use StringIO to convert the string output to a file-like object
                blast_output_io = io.StringIO(blast_output)
                blast_records = NCBIXML.parse(blast_output_io)
                for blast_record in blast_records:
                    for alignment in blast_record.alignments:
                        for hsp in alignment.hsps:
                            evalue = hsp.expect
                            evalues.append(evalue)
                            results += f"Query Sequence: {seq_record.id}\n"
                            results += f"Alignment Title: {alignment.title}\n"
                            results += f"Length: {alignment.length}\n"
                            results += f"E-value: {evalue}\n"
                            results += f"Score: {hsp.score}\n\n"

                            # Check if Cry gene is detected
                            if any(cry_gene_file in alignment.title for cry_gene_file in fasta_files):
                                cry_gene_found = True

                        # Provide specific message based on Cry gene detection
                        if cry_gene_found:
                            results += f"Plant is likely a GMO (contains Cry gene)\n"
                        else:
                            results += f"Plant doesn't contain Cry gene\n"

            result_handler.display_results(results, evalues)

        blast_thread = threading.Thread(target=process_blast_results,
                                        args=(query_sequences, selected_cry_genes_path, evalue_threshold, word_size, gap_open, gap_extend, result_handler))
        blast_thread.start()

    except Exception as e:
        result_handler.handle_error(str(e))
