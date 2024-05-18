# blastn.py
import os
import wx
import threading
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

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

            for file_name in os.listdir(cry_genes_path):
                if file_name.endswith(".fasta"):
                    cry_gene_path = os.path.join(cry_genes_path, file_name)
                    cry_gene_name = os.path.basename(cry_gene_path).replace(".fasta", "")
                    cry_gene_content = str(SeqIO.read(cry_gene_path, "fasta").seq)

                    for seq_record in query_sequences:
                        query_gene = str(seq_record.seq)

                        result_handle = NCBIWWW.qblast("blastn", "nr", query_gene, expect=evalue_threshold, word_size=word_size, gapcosts=(gap_open, gap_extend))
                        blast_records = NCBIXML.parse(result_handle)

                        cry_gene_found = False

                        for record in blast_records:
                            for alignment in record.alignments:
                                evalue = alignment.hsps[0].expect
                                evalues.append(evalue)
                                results += f"Query Sequence: {seq_record.id}\n"
                                results += f"Alignment Title: {alignment.title}\n"
                                results += f"Length: {alignment.length}\n"
                                results += f"E-value: {evalue}\n"
                                results += f"Score: {alignment.hsps[0].score}\n\n"

                                # Check if Cry gene is detected
                                if cry_gene_content in alignment.title:
                                    cry_gene_found = True

                        # Provide specific message based on Cry gene detection
                        if cry_gene_found:
                            results += f"Plant is likely a GMO (contains {cry_gene_name})\n"
                        else:
                            results += f"Plant doesn't contain {cry_gene_name}\n"

            result_handler.display_results(results, evalues)

        blast_thread = threading.Thread(target=process_blast_results,
                                        args=(query_sequences, selected_cry_genes_path, evalue_threshold, word_size, gap_open, gap_extend, result_handler))
        blast_thread.start()

    except Exception as e:
        result_handler.handle_error(str(e))
